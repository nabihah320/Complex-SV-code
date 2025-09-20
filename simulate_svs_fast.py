#!/usr/bin/env python3
"""
Fast SV Simulator with improved repeat anchoring
- Scans repeats across the whole chromosome quickly (sampled windows, merged)
- Supports anchoring SVs to repeats by start/span/breakpoints
- Guarantees BND and INV are NOT in repeat regions
"""

import random
import sys
import numpy as np
from collections import defaultdict, OrderedDict

# ---------------------------
# Utilities
# ---------------------------

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, b) for b in seq[::-1])

def merge_intervals(intervals, pad=0):
    """Merge intervals (start, end, type). Keeps the predominant type when merging."""
    if not intervals:
        return []
    # Sort by start
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    cur_s, cur_e, cur_t = intervals[0]
    for s, e, t in intervals[1:]:
        if s <= cur_e + pad:
            # merge
            new_s = cur_s
            new_e = max(cur_e, e)
            # keep type if same, otherwise mark generic "repeat"
            new_t = cur_t if cur_t == t else "repeat"
            cur_s, cur_e, cur_t = new_s, new_e, new_t
        else:
            merged.append((cur_s, cur_e, cur_t))
            cur_s, cur_e, cur_t = s, e, t
    merged.append((cur_s, cur_e, cur_t))
    return merged

def interval_overlaps(intervals, start, end):
    """Return True if [start, end) overlaps any interval; also if start or end inside."""
    # intervals must be sorted/merged
    lo, hi = 0, len(intervals) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e, _ = intervals[mid]
        if end <= s:
            hi = mid - 1
        elif start >= e:
            lo = mid + 1
        else:
            return True
    return False

def point_in_interval(intervals, pos):
    """Return True if pos is inside any interval."""
    lo, hi = 0, len(intervals) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e, _ = intervals[mid]
        if pos < s:
            hi = mid - 1
        elif pos >= e:
            lo = mid + 1
        else:
            return True
    return False

# ---------------------------
# Main
# ---------------------------

def simulate_svs_fast(ref_file, output_dir):
    random.seed(42)
    np.random.seed(42)

    # Load reference
    print("Loading reference...")
    sequences = OrderedDict()
    chrom_lengths = OrderedDict()
    current_seq = []
    current_chr = None

    with open(ref_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_chr:
                    seq = ''.join(current_seq)
                    sequences[current_chr] = seq
                    chrom_lengths[current_chr] = len(seq)
                current_chr = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        if current_chr:
            seq = ''.join(current_seq)
            sequences[current_chr] = seq
            chrom_lengths[current_chr] = len(seq)

    print(f"Loaded {len(sequences)} chromosomes")

    # Fast repeat detection across ENTIRE chromosome (sampled)
    print("Fast repeat region detection...")
    sampled_step = 10000     # sample every 10 kb
    window_len   = 600       # 600 bp window
    repeat_regions_raw = defaultdict(list)

    for chrom, seq in sequences.items():
        seq_len = len(seq)
        for i in range(0, max(0, seq_len - window_len), sampled_step):
            window = seq[i:i+window_len]

            # Homopolymer
            is_repeat = False
            for base in ['A', 'T', 'G', 'C']:
                if base * 12 in window:
                    repeat_regions_raw[chrom].append((i, i+window_len, 'homopolymer'))
                    is_repeat = True
                    break

            if is_repeat:
                continue

            # Quick tandem checks (2-3bp)
            for size in [2, 3]:
                # look for unit repeated 4+ times anywhere near the start
                unit = window[:size]
                if unit * 4 == window[:size*4]:
                    repeat_regions_raw[chrom].append((i, i+window_len, 'tandem'))
                    is_repeat = True
                    break

    # Merge nearby repeat windows for better coverage
    repeat_regions = {}
    total_repeats = 0
    for chrom, ivals in repeat_regions_raw.items():
        merged = merge_intervals(ivals, pad=500)  # merge windows with <=500bp gap
        repeat_regions[chrom] = merged
        total_repeats += len(merged)
    print(f"  Found {total_repeats} merged repeat regions (fast scan)")

    # Helper: check repeat relationship according to anchor policy
    def in_repeat_by_policy(chrom, pos, size, anchor='span'):
        ivals = repeat_regions.get(chrom, [])
        if not ivals:
            return False
        start = pos
        end = pos + max(1, size)
        if anchor == 'start':
            return point_in_interval(ivals, start)
        elif anchor == 'breakpoints':
            return point_in_interval(ivals, start) or point_in_interval(ivals, end-1)
        else:  # 'span'
            return interval_overlaps(ivals, start, end)

    # Conflict tracking
    sv_positions = defaultdict(list)  # per chrom list of (start, end)

    def conflict_exists(chrom, start, end, pad=500):
        for s, e in sv_positions[chrom]:
            if not (end < s - pad or start > e + pad):
                return True
        return False

    # Position finder
    def find_position(chrom, size, require_repeat=False, forbid_repeat=False, anchor='span',
                      max_attempts=120, targeted_attempts=30):
        """
        If require_repeat: try targeting repeat-anchored placement first, then random but must satisfy repeat policy.
        If forbid_repeat: ensure NOT in repeat by the given anchor policy.
        anchor: 'start', 'span', or 'breakpoints'
        """
        chrom_len = len(sequences[chrom])
        ivals = repeat_regions.get(chrom, [])

        # Targeted repeat anchoring
        if require_repeat and ivals:
            for _ in range(targeted_attempts):
                rs, re, _ = random.choice(ivals)
                # Choose pos so anchor policy is satisfied without requiring whole span inside repeat
                if anchor == 'start':
                    pos = random.randint(max(50, rs), max(50, min(re-1, chrom_len-2)))
                    pos = min(pos, chrom_len - size - 50)
                elif anchor == 'breakpoints':
                    # Place so either start or end falls inside the repeat
                    # Try to set start inside repeat first
                    pos = random.randint(max(50, rs), max(50, min(re-1, chrom_len-2)))
                    pos = min(pos, chrom_len - size - 50)
                    if not in_repeat_by_policy(chrom, pos, size, 'breakpoints'):
                        # Try to set end inside repeat: start somewhere before repeat
                        pos = random.randint(max(50, rs - size + 1), max(50, rs))
                        pos = max(50, min(pos, chrom_len - size - 50))
                else:  # 'span' -> any overlap
                    pos = random.randint(max(50, rs - size + 1), max(50, min(re-1, chrom_len-2)))
                    pos = max(50, min(pos, chrom_len - size - 50))

                if pos < 50 or pos + size > chrom_len - 50:
                    continue
                if conflict_exists(chrom, pos, pos + size):
                    continue
                if in_repeat_by_policy(chrom, pos, size, anchor):
                    sv_positions[chrom].append((pos, pos + size))
                    return pos, True

        # Random attempts
        for _ in range(max_attempts):
            pos = random.randint(5000, max(5001, chrom_len - size - 5000))
            if conflict_exists(chrom, pos, pos + size):
                continue
            in_rep = in_repeat_by_policy(chrom, pos, size, anchor)
            if require_repeat and not in_rep:
                continue
            if forbid_repeat and in_rep:
                continue
            sv_positions[chrom].append((pos, pos + size))
            return pos, in_rep

        # Final relaxed attempts if require_repeat fails: accept any, but mark repeat status
        if require_repeat:
            for _ in range(40):
                pos = random.randint(5000, max(5001, chrom_len - size - 5000))
                if conflict_exists(chrom, pos, pos + size):
                    continue
                in_rep = in_repeat_by_policy(chrom, pos, size, anchor)
                sv_positions[chrom].append((pos, pos + size))
                return pos, in_rep

        return None, False

    # ---------------------------
    # Generate SVs
    # ---------------------------

    svs = []
    sv_id = 1
    chroms = list(sequences.keys())

    sv_counts = {
        'DEL': 150,
        'INS': 120,
        'DUP': 60,
        'INV': 40,
        'BND': 30,
        'COMPLEX': 20
    }

    # DEL
    print(f"Generating {sv_counts['DEL']} deletions...")
    for i in range(sv_counts['DEL']):
        chrom = random.choice(chroms)
        size = int(np.random.lognormal(6.0, 1.5))
        size = max(50, min(100000, size))
        prefer_repeat = (i % 3 == 0)  # ~33% try to anchor to repeats
        pos, in_repeat = find_position(
            chrom, size,
            require_repeat=prefer_repeat,
            forbid_repeat=False,
            anchor='span'
        )
        if pos is not None:
            support = random.randint(20, 40) if in_repeat else random.randint(15, 35)
            qual = random.randint(25, 50) if in_repeat else random.randint(30, 60)
            svs.append({
                'chrom': chrom, 'pos': pos, 'end': pos + size, 'type': 'DEL',
                'size': size, 'id': f'DEL_{sv_id}',
                'genotype': '0/1' if random.random() < 0.7 else '1/1',
                'support': support, 'af': 0.5 if random.random() < 0.7 else 1.0,
                'qual': qual, 'in_repeat': in_repeat
            })
            sv_id += 1

    # INS
    print(f"Generating {sv_counts['INS']} insertions...")
    for i in range(sv_counts['INS']):
        chrom = random.choice(chroms)
        size = int(np.random.lognormal(5.0, 1.2))
        size = max(50, min(5000, size))
        prefer_repeat = (i % 2 == 0)  # ~50% in repeats
        pos, in_repeat = find_position(
            chrom, max(1, size // 10),  # tiny span for placement; anchor at start
            require_repeat=prefer_repeat,
            forbid_repeat=False,
            anchor='start'
        )
        if pos is not None:
            if in_repeat and random.random() < 0.7:
                repeat_unit = random.choice(['A', 'T', 'AT', 'GC', 'CAG'])
                seq_ins = (repeat_unit * (size // len(repeat_unit) + 2))[:size]
            else:
                seq_ins = ''.join(random.choices('ACGT', k=size))
            support = random.randint(15, 35) if in_repeat else random.randint(12, 30)
            qual = random.randint(20, 45) if in_repeat else random.randint(25, 55)
            svs.append({
                'chrom': chrom, 'pos': pos, 'end': pos, 'type': 'INS',
                'size': size, 'id': f'INS_{sv_id}', 'seq': seq_ins,
                'genotype': '0/1' if random.random() < 0.8 else '1/1',
                'support': support, 'af': 0.5 if random.random() < 0.8 else 1.0,
                'qual': qual, 'in_repeat': in_repeat
            })
            sv_id += 1

    # DUP (tandem + interspersed)
    print(f"Generating {sv_counts['DUP']} duplications...")
    tandem_count = int(sv_counts['DUP'] * 0.7)
    for i in range(sv_counts['DUP']):
        chrom = random.choice(chroms)
        size = int(np.random.lognormal(7.0, 1.5))
        size = max(100, min(20000, size))
        is_tandem = i < tandem_count

        # For tandems, anchor breakpoints to repeats ~60% of the time
        want_repeat = is_tandem and (random.random() < 0.6)
        pos, in_repeat = find_position(
            chrom, size,
            require_repeat=want_repeat,
            forbid_repeat=False,
            anchor='breakpoints'
        )
        if pos is not None:
            support = random.randint(12, 28) if in_repeat else random.randint(10, 25)
            qual = random.randint(20, 45) if in_repeat else random.randint(20, 50)
            sv_data = {
                'chrom': chrom, 'pos': pos, 'end': pos + size, 'type': 'DUP',
                'subtype': 'TANDEM' if is_tandem else 'INTERSPERSED',
                'size': size, 'id': f'DUP_{sv_id}',
                'genotype': '0/1', 'cn': 2,
                'support': support, 'af': 0.5,
                'qual': qual, 'in_repeat': in_repeat
            }
            if not is_tandem:
                # Interspersed target position on same chrom, away from source
                chrom_len = len(sequences[chrom])
                # pick a target at least 50 kb away, keep in bounds
                offset = random.randint(50000, 200000)
                target = pos + offset
                if target >= chrom_len - 1000:
                    target = max(1000, pos - offset)
                sv_data['target_pos'] = max(1000, min(chrom_len - 1000, target))
            svs.append(sv_data)
            sv_id += 1

    # INV — ensure NOT in repeats
    print(f"Generating {sv_counts['INV']} inversions...")
    for i in range(sv_counts['INV']):
        chrom = random.choice(chroms)
        size = int(np.random.lognormal(7.5, 1.8))
        size = max(200, min(50000, size))
        pos, in_repeat = find_position(
            chrom, size,
            require_repeat=False,
            forbid_repeat=True,     # ensure NOT in repeat
            anchor='span'
        )
        if pos is not None:
            support = random.randint(10, 25)
            qual = random.randint(25, 55)
            svs.append({
                'chrom': chrom, 'pos': pos, 'end': pos + size, 'type': 'INV',
                'size': size, 'id': f'INV_{sv_id}',
                'genotype': '0/1' if random.random() < 0.9 else '1/1',
                'support': support, 'af': 0.5 if random.random() < 0.9 else 1.0,
                'qual': qual, 'in_repeat': False   # guaranteed by forbid_repeat
            })
            sv_id += 1

    # BND — ensure NOT in repeats
    print(f"Generating {sv_counts['BND']} breakends...")
    for i in range(sv_counts['BND']):
        if len(chroms) < 1:
            break
        # choose chrom1; ~70% inter-chrom as before, but both ends must avoid repeats
        chrom1 = random.choice(chroms)
        pos1, in_rep1 = find_position(
            chrom1, 1,
            require_repeat=False,
            forbid_repeat=True,     # not in repeat
            anchor='start'
        )
        if pos1 is None:
            continue

        inter_chrom = len(chroms) >= 2 and random.random() < 0.7
        if inter_chrom:
            chrom2 = random.choice([c for c in chroms if c != chrom1])
            # find position not in repeats on chrom2
            pos2, in_rep2 = find_position(
                chrom2, 1,
                require_repeat=False,
                forbid_repeat=True,  # not in repeat
                anchor='start'
            )
            if pos2 is None:
                continue
        else:
            chrom2 = chrom1
            # different position on same chrom, not in repeat
            pos2, in_rep2 = find_position(
                chrom2, 1,
                require_repeat=False,
                forbid_repeat=True,
                anchor='start'
            )
            if pos2 is None:
                continue

        support = random.randint(5, 15)
        qual = random.randint(15, 40)
        orientation = random.choice(['N[', 'N]', '[N', ']N'])
        svs.append({
            'chrom': chrom1, 'pos': pos1, 'end': pos1, 'type': 'BND',
            'size': 0, 'id': f'BND_{sv_id}',
            'chr2': chrom2, 'pos2': pos2, 'orientation': orientation,
            'genotype': '0/1', 'support': support, 'af': random.uniform(0.2, 0.5),
            'qual': qual, 'in_repeat': False
        })
        sv_id += 1

    # CPX — anchor to repeats for ~50% using breakpoint policy
    print(f"Generating {sv_counts['COMPLEX']} complex SVs...")
    patterns = [['DEL', 'INV'], ['DUP', 'INV'], ['DEL', 'INS'], ['DEL', 'DUP']]
    for i in range(sv_counts['COMPLEX']):
        chrom = random.choice(chroms)
        pattern = random.choice(patterns)
        total_size = random.randint(1000, 8000)
        want_repeat = (random.random() < 0.5)
        pos, in_repeat = find_position(
            chrom, total_size,
            require_repeat=want_repeat,
            forbid_repeat=False,
            anchor='breakpoints'
        )
        if pos is not None:
            support = random.randint(10, 25) if in_repeat else random.randint(8, 20)
            qual = random.randint(18, 40) if in_repeat else random.randint(20, 45)
            # Build components left-to-right; INV inside CPX may cross repeats; we keep the CPX anchor policy only for the outer CPX
            current_pos = pos
            component_size = max(100, total_size // len(pattern))
            for comp_type in pattern:
                comp_start = current_pos
                comp_end = current_pos + component_size
                comp_sv = {
                    'chrom': chrom, 'pos': comp_start, 'end': comp_end,
                    'type': comp_type, 'size': component_size,
                    'id': f'COMPLEX_{sv_id}_part{comp_type}',
                    'complex_id': f'COMPLEX_{sv_id}',
                    'genotype': '0/1',
                    'support': support, 'af': 0.5,
                    'qual': qual, 'in_repeat': in_repeat
                }
                if comp_type == 'INS':
                    comp_sv['seq'] = ''.join(random.choices('ACGT', k=component_size))
                    comp_sv['end'] = comp_start
                elif comp_type == 'DUP':
                    comp_sv['subtype'] = 'TANDEM'
                    comp_sv['cn'] = 2
                # Note: If you want inner INV strictly non-repeat, you could enforce here,
                # but the user's requirement only forbids top-level INV SVs, not CPX parts.
                svs.append(comp_sv)
                current_pos = comp_end + random.randint(50, 200)

            svs.append({
                'chrom': chrom, 'pos': pos, 'end': pos + total_size, 'type': 'CPX',
                'size': total_size, 'id': f'COMPLEX_{sv_id}', 'components': ','.join(pattern),
                'genotype': '0/1', 'support': support, 'af': 0.5,
                'qual': qual, 'in_repeat': in_repeat
            })
            sv_id += 1

    # ---------------------------
    # Write VCF
    # ---------------------------

    print("Writing truth VCF...")
    with open(f"{output_dir}/truth.vcf", 'w') as f:
        # Headers
        f.write("##fileformat=VCFv4.3\n")
        f.write(f"##reference={ref_file}\n")
        for chrom, length in chrom_lengths.items():
            f.write(f"##contig=<ID={chrom},length={length}>\n")

        # ALT headers
        f.write('##ALT=<ID=DEL,Description="Deletion">\n')
        f.write('##ALT=<ID=INS,Description="Insertion">\n')
        f.write('##ALT=<ID=DUP,Description="Duplication">\n')
        f.write('##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">\n')
        f.write('##ALT=<ID=DUP:INTERSPERSED,Description="Interspersed Duplication">\n')
        f.write('##ALT=<ID=INV,Description="Inversion">\n')
        f.write('##ALT=<ID=BND,Description="Breakend">\n')
        f.write('##ALT=<ID=CPX,Description="Complex SV">\n')

        # INFO fields
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        f.write('##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Read support">\n')
        f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n')
        f.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome 2 for BND">\n')
        f.write('##INFO=<ID=CN,Number=1,Type=Float,Description="Copy number">\n')
        f.write('##INFO=<ID=DUPTYPE,Number=1,Type=String,Description="Duplication type">\n')
        f.write('##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Complex SV components">\n')
        f.write('##INFO=<ID=REPEAT,Number=1,Type=String,Description="In repeat region">\n')

        # FORMAT
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        # Write variants
        for sv in sorted(svs, key=lambda x: (x['chrom'], x['pos'])):
            chrom = sv['chrom']
            pos = sv['pos'] + 1
            sv_id_str = sv['id']
            qual = sv.get('qual', 30)

            ref = 'N'
            if sv['type'] == 'DEL':
                alt = '<DEL>'
                info = f"SVTYPE=DEL;END={sv['end']+1};SVLEN=-{sv['size']};SUPPORT={sv['support']};AF={sv['af']:.3f}"
            elif sv['type'] == 'INS':
                alt = '<INS>'
                info = f"SVTYPE=INS;END={sv['end']+1};SVLEN={sv['size']};SUPPORT={sv['support']};AF={sv['af']:.3f}"
            elif sv['type'] == 'DUP':
                subtype = sv.get('subtype', 'TANDEM')
                alt = f"<DUP:{subtype}>"
                info = f"SVTYPE=DUP;END={sv['end']+1};SVLEN={sv['size']};SUPPORT={sv['support']};AF={sv['af']:.3f};DUPTYPE={subtype};CN={sv.get('cn', 2)}"
            elif sv['type'] == 'INV':
                alt = '<INV>'
                info = f"SVTYPE=INV;END={sv['end']+1};SVLEN={sv['size']};SUPPORT={sv['support']};AF={sv['af']:.3f}"
            elif sv['type'] == 'BND':
                chr2 = sv.get('chr2', chrom)
                pos2 = sv.get('pos2', pos) + 1
                orientation = sv.get('orientation', 'N[')
                if orientation == 'N[':
                    alt = f"N[{chr2}:{pos2}["
                elif orientation == 'N]':
                    alt = f"N]{chr2}:{pos2}]"
                elif orientation == '[N':
                    alt = f"[{chr2}:{pos2}[N"
                else:
                    alt = f"]{chr2}:{pos2}]N"
                info = f"SVTYPE=BND;SUPPORT={sv['support']};AF={sv['af']:.3f};CHR2={chr2}"
            elif sv['type'] == 'CPX':
                alt = '<CPX>'
                info = f"SVTYPE=CPX;END={sv['end']+1};SVLEN={sv['size']};SUPPORT={sv['support']};AF={sv['af']:.3f};CPX_TYPE={sv.get('components', 'Unknown')}"
            else:
                continue

            if sv.get('in_repeat', False):
                info += ";REPEAT=Yes"
            else:
                info += ";REPEAT=No"

            genotype = sv.get('genotype', '0/1')
            f.write(f"{chrom}\t{pos}\t{sv_id_str}\t{ref}\t{alt}\t{qual}\tPASS\t{info}\tGT\t{genotype}\n")

    # ---------------------------
    # Apply SVs to reference
    # ---------------------------

    print("Applying SVs to reference...")
    modified_seqs = {}
    for chrom, seq in sequences.items():
        seq_list = list(seq)
        chrom_svs = sorted([sv for sv in svs if sv['chrom'] == chrom and sv['type'] not in ['BND', 'CPX']],
                           key=lambda x: x['pos'], reverse=True)
        for sv in chrom_svs:
            try:
                if sv['type'] == 'DEL':
                    del seq_list[sv['pos']:sv['end']]
                elif sv['type'] == 'INS':
                    seq_list[sv['pos']:sv['pos']] = list(sv.get('seq', 'N' * sv['size']))
                elif sv['type'] == 'DUP':
                    dup_seq = seq_list[sv['pos']:sv['end']]
                    if sv.get('subtype') == 'TANDEM':
                        seq_list[sv['end']:sv['end']] = dup_seq
                    else:
                        target = sv.get('target_pos', sv['end'] + 50000)
                        if 0 < target < len(seq_list):
                            seq_list[target:target] = dup_seq
                elif sv['type'] == 'INV':
                    inv_seq = ''.join(seq_list[sv['pos']:sv['end']])
                    seq_list[sv['pos']:sv['end']] = list(reverse_complement(inv_seq))
            except Exception:
                pass
        modified_seqs[chrom] = ''.join(seq_list)

    with open(f"{output_dir}/simulated.fasta", 'w') as f:
        for chrom, seq in modified_seqs.items():
            f.write(f">{chrom}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    # ---------------------------
    # Summary
    # ---------------------------

    sv_type_counts = defaultdict(int)
    repeat_counts = defaultdict(int)
    for sv in svs:
        sv_type_counts[sv['type']] += 1
        if sv.get('in_repeat', False):
            repeat_counts[sv['type']] += 1

    print(f"\nGenerated {len(svs)} SVs:")
    for sv_type in ['DEL', 'INS', 'DUP', 'INV', 'BND', 'CPX']:
        count = sv_type_counts[sv_type]
        repeat_count = repeat_counts.get(sv_type, 0)
        if count > 0:
            print(f"  {sv_type}: {count} ({repeat_count} in repeats)")

if __name__ == "__main__":
    simulate_svs_fast(sys.argv[1], sys.argv[2])
