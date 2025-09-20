#!/usr/bin/env python3
"""
Complex_SV Caller with Sniffles2-inspired filtering and repeat handling
Version 3.1 - Complete Implementation with all fixes
"""

import pysam
import sys
import json
import numpy as np
from collections import defaultdict, Counter, namedtuple
from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional, Tuple
import logging
from scipy import stats
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

Evidence = namedtuple('Evidence', ['type', 'chrom', 'start', 'end', 'sv_type', 'size', 'read_name', 'mapq', 'extra'])

@dataclass
class SVSignature:
    chrom: str
    start: int
    end: int
    sv_type: str
    size: int
    support: int = 1
    mapq: int = 20
    reads: Set[str] = field(default_factory=set)
    source: str = 'unknown'
    af: float = 0.0
    cn: float = 1.0
    std_pos: float = 0.0
    std_len: float = 0.0
    precise: bool = False
    re: int = 0
    rnames: List[str] = field(default_factory=list)
    mate_chr: Optional[str] = None
    mate_pos: Optional[int] = None
    bnd_type: Optional[str] = None
    split_support: int = 0
    spanning_support: int = 0
    depth_support: int = 0
    genotype: str = './.'
    genotype_qual: int = 0
    ref_support: int = 0
    alt_support: int = 0
    is_complex: bool = False
    complex_type: Optional[str] = None
    complex_components: List[str] = field(default_factory=list)
    filter: str = 'PASS'
    qual: float = 0.0
    strand_support: Dict[str, int] = field(default_factory=dict)
    is_artifact: bool = False
    in_repeat: bool = False
    repeat_type: Optional[str] = None
    ref_context: Optional[str] = None

class RepeatHandler:
    """Handle repeat regions and their impact on SV calling"""
    
    def __init__(self, ref_file=None, trf_bed=None):
        self.ref = pysam.FastaFile(ref_file) if ref_file else None
        self.trf_regions = self.load_trf_bed(trf_bed) if trf_bed else {}
        
        # Repeat detection parameters
        self.homopolymer_threshold = 7
        self.tandem_threshold = 0.8
        self.low_complexity_threshold = 0.3
        
    def load_trf_bed(self, trf_bed):
        """Load tandem repeat regions from TRF bed file"""
        trf_regions = defaultdict(list)
        try:
            with open(trf_bed, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        repeat_info = parts[3] if len(parts) > 3 else 'tandem'
                        trf_regions[chrom].append((start, end, repeat_info))
        except:
            logger.warning(f"Could not load TRF bed file: {trf_bed}")
        return trf_regions
    
    def check_repeat_context(self, chrom, pos, end):
        """Check if position is in a repeat region"""
        repeat_info = {
            'in_repeat': False,
            'repeat_type': None,
            'homopolymer': False,
            'tandem': False,
            'low_complexity': False
        }
        
        # First, check TRF regions (most reliable)
        if chrom in self.trf_regions:
            for trf_start, trf_end, trf_type in self.trf_regions[chrom]:
                if not (end < trf_start or pos > trf_end):
                    repeat_info['in_repeat'] = True
                    repeat_info['tandem'] = True
                    repeat_info['repeat_type'] = 'tandem'
                    return repeat_info  # Trust TRF, return immediately
        
        # Only check sequence if no TRF annotation
        if self.ref and chrom in self.ref.references:
            try:
                context_start = max(0, pos - 20)  # Smaller window
                context_end = min(self.ref.get_reference_length(chrom), end + 20)
                context = self.ref.fetch(chrom, context_start, context_end).upper()
                
                # Stricter homopolymer detection
                for base in ['A', 'T', 'G', 'C']:
                    if base * 12 in context:  # Increased from 8 to 12
                        repeat_info['homopolymer'] = True
                        repeat_info['in_repeat'] = True
                        repeat_info['repeat_type'] = 'homopolymer'
                        break
                

            except:
                pass
        
        return repeat_info
    
    def adjust_thresholds_for_repeat(self, sv_call, base_threshold):
        """Adjust filtering thresholds based on repeat context"""
        if sv_call.in_repeat:
            if sv_call.repeat_type == 'homopolymer':
                # Need more support in homopolymers
                if sv_call.sv_type in ['INS', 'DEL']:
                    return base_threshold * 1.5
            elif sv_call.repeat_type == 'tandem':
                # Tandem repeats prone to artifacts
                if sv_call.sv_type == 'DUP':
                    return base_threshold * 0.8  # More lenient for expected DUPs
                else:
                    return base_threshold * 1.3
            elif sv_call.repeat_type == 'low_complexity':
                return base_threshold * 1.2
        return base_threshold

class SVQualityScorer:
    """Sniffles2-inspired quality scoring and filtering"""
    
    def __init__(self, coverage, ploidy=2):
        self.coverage = coverage
        self.ploidy = ploidy
        self.expected_alt_support = self.coverage / self.ploidy
        
        # Sniffles2-style thresholds
        self.min_qual_threshold = 10
        self.strand_bias_threshold = 0.1
        self.position_stdev_threshold = 0.5
        
    def calculate_sv_quality(self, sv_call):
        """Calculate quality score similar to Sniffles2"""
        
        # Base quality from support
        base_qual = self.calculate_support_quality(sv_call)
        
        # Adjust for strand bias
        strand_qual = self.calculate_strand_quality(sv_call)
        
        # Adjust for position consistency
        position_qual = self.calculate_position_quality(sv_call)
        
        # Adjust for AF consistency
        af_qual = self.calculate_af_quality(sv_call)
        
        # Combine scores (Sniffles2 approach)
        final_qual = (base_qual * 0.4 + 
                     strand_qual * 0.2 + 
                     position_qual * 0.2 + 
                     af_qual * 0.2)
        
        # Apply SV-type specific adjustments
        if sv_call.sv_type == 'BND':
            if sv_call.split_support > 0:
                final_qual *= 1.2
        elif sv_call.sv_type == 'DUP':
            if sv_call.depth_support > 0:
                final_qual *= 1.3
        elif sv_call.sv_type == 'INS':
            if sv_call.std_len < 10:
                final_qual *= 1.1
        
        # Penalty for repeat regions
        if sv_call.in_repeat:
            if sv_call.repeat_type == 'homopolymer':
                final_qual *= 0.7
            elif sv_call.repeat_type == 'low_complexity':
                final_qual *= 0.8
        
        return min(99, final_qual)
    
    def calculate_support_quality(self, sv_call):
        """Calculate quality based on read support (Sniffles2 method)"""
        if sv_call.support == 0:
            return 0
        
        # Expected support for heterozygous variant
        if sv_call.genotype == '0/1':
            expected = self.expected_alt_support
        elif sv_call.genotype == '1/1':
            expected = self.expected_alt_support * 2
        else:
            expected = self.expected_alt_support
        
        # Poisson probability
        if expected > 0:
            if sv_call.support >= expected:
                qual = -10 * np.log10(1 - stats.poisson.cdf(sv_call.support - 1, expected) + 1e-10)
            else:
                qual = -10 * np.log10(stats.poisson.cdf(sv_call.support, expected) + 1e-10)
        else:
            qual = sv_call.support * 3
        
        return min(60, qual)
    
    def calculate_strand_quality(self, sv_call):
        """Calculate strand bias penalty (Sniffles2 approach)"""
        if not sv_call.strand_support:
            return 30
        
        forward = sv_call.strand_support.get('forward', 0)
        reverse = sv_call.strand_support.get('reverse', 0)
        total = forward + reverse
        
        if total == 0:
            return 0
        
        strand_ratio = min(forward, reverse) / total
        
        if strand_ratio < self.strand_bias_threshold:
            return 10
        else:
            return 30
    
    def calculate_position_quality(self, sv_call):
        """Calculate position consistency quality"""
        if sv_call.std_pos == 0:
            return 40
        
        relative_stdev = sv_call.std_pos / max(1, sv_call.size)
        
        if relative_stdev < 0.1:
            return 35
        elif relative_stdev < 0.3:
            return 25
        elif relative_stdev < 0.5:
            return 15
        else:
            return 5
    
    def calculate_af_quality(self, sv_call):
        """Calculate allele frequency quality"""
        if sv_call.af < 0.1:
            return 5
        elif sv_call.af < 0.2:
            return 15
        elif sv_call.af < 0.35:
            return 25
        else:
            return 30

class ArtifactFilter:
    """Sniffles2-inspired artifact filtering"""
    
    def __init__(self, repeat_handler=None):
        self.repeat_handler = repeat_handler
        
    def is_likely_artifact(self, sv_call, reads_in_region):
        """Check if SV is likely an artifact (Sniffles2 approach)"""
        
        # Don't mark DUPs as artifacts if they have depth evidence
        if sv_call.sv_type == 'DUP' and (sv_call.source == 'depth' or sv_call.depth_support > 0):
            return False
        
        # Check repeat context
        if self.repeat_handler:
            repeat_info = self.repeat_handler.check_repeat_context(
                sv_call.chrom, sv_call.start, sv_call.end
            )
            
            # In homopolymer with low support = likely artifact
            if repeat_info['homopolymer'] and sv_call.sv_type in ['INS']:
                if sv_call.support < 5:
                    return True
            
            if repeat_info['homopolymer'] and sv_call.sv_type in ['DEL']:
                if sv_call.support < 1:
                    return True

            # Low complexity with low AF = likely artifact
            if repeat_info['low_complexity'] and sv_call.af < 0.2:
                return True
        
        # Check for alignment artifacts
        if self.has_alignment_artifact_signature(sv_call, reads_in_region):
            return True
        
        # Check for strand bias artifacts
        if self.has_extreme_strand_bias(sv_call):
            return True
        
        return False
    
    def has_alignment_artifact_signature(self, sv_call, reads_in_region):
        """Detect alignment artifacts (Sniffles2 approach)"""
        if not reads_in_region:
            return False
        
        # BNDs with all reads having same low MAPQ might be artifacts
        if sv_call.sv_type == 'BND':
            mapqs = [r.mapping_quality for r in reads_in_region if not r.is_supplementary]
            if mapqs and len(set(mapqs)) == 1 and mapqs[0] < 20:
                return True
        
        # Too many supplementary alignments
        if reads_in_region:
            supp_count = sum(1 for r in reads_in_region if r.is_supplementary)
            if supp_count > len(reads_in_region) * 0.8 and sv_call.support < 5:
                return True
        
        return False
    
    def has_extreme_strand_bias(self, sv_call):
        """Check for extreme strand bias indicating artifact"""
        if sv_call.strand_support:
            forward = sv_call.strand_support.get('forward', 0)
            reverse = sv_call.strand_support.get('reverse', 0)
            total = forward + reverse
            
            if total >= 4:
                ratio = min(forward, reverse) / total
                if ratio < 0.05:
                    return True
        return False

class EnhancedSVCaller:
    def __init__(self, bam_file, ref_file, params_file, trf_bed=None):
        self.bam = pysam.AlignmentFile(bam_file, "rb")
        self.ref = pysam.FastaFile(ref_file) if ref_file else None
        
        with open(params_file) as f:
            self.params = json.load(f)
        
        self.min_sv_size = self.params.get('min_sv_size', 50)
        self.min_mapq = self.params.get('min_mapq', 10)
        self.min_af = self.params.get('min_af', 0.10)
        self.expected_coverage = float(self.params.get('expected_coverage', 30))
        
        # Complex SV parameters
        self.complex_min_signatures = self.params.get('complex_min_signatures', 2)
        self.complex_max_distance = self.params.get('complex_max_distance', 1000)
        self.complex_merge_distance = self.params.get('complex_merge_distance', 500)
        
        # Clustering parameters
        self.cluster_max_distance = 0.3
        self.min_support = 2
        self.min_support_bnd = self.params.get('min_support_bnd', 1)
        self.min_support_complex = self.params.get('min_support_complex', 2)
        
        # Store reference info
        self.ref_lengths = {}
        for i in range(self.bam.nreferences):
            self.ref_lengths[self.bam.references[i]] = self.bam.lengths[i]
        
        # Initialize repeat handler
        self.repeat_handler = RepeatHandler(ref_file, trf_bed)
        
        # Estimate coverage
        self.coverage = self.estimate_coverage_properly()
        self.ploidy = self.params.get('ploidy', 2)
        logger.info(f"Estimated coverage: {self.coverage:.1f}x")
        
        # Initialize scorers and filters
        self.quality_scorer = SVQualityScorer(self.coverage, self.ploidy)
        self.artifact_filter = ArtifactFilter(self.repeat_handler)
        
        # Evidence collection
        self.evidence = defaultdict(list)
        self.segments = defaultdict(list)
        self.bnd_evidence = []
        self.interspersed_evidence = []
        self.complex_signatures = []
        
    def estimate_coverage_properly(self):
        """More robust coverage estimation"""
        coverage_samples = []
        
        logger.info("Estimating coverage from BAM file...")
        
        # Try multiple approaches
        total_reads = 0
        total_bases = 0
        
        # Method 1: Sample multiple chromosomes and positions
        for ref_idx in range(min(5, self.bam.nreferences)):
            chrom = self.bam.references[ref_idx]
            chrom_len = self.bam.lengths[ref_idx]
            
            if chrom_len < 10000:
                continue
            
            # Sample 20 regions of 1kb each
            for sample_num in range(20):
                try:
                    # Sample from different parts of chromosome
                    start_pos = sample_num * (chrom_len // 20) + 1000
                    end_pos = min(start_pos + 1000, chrom_len - 1000)
                    
                    if end_pos <= start_pos:
                        continue
                    
                    region_depth = 0
                    region_length = 0
                    
                    for pileup_col in self.bam.pileup(chrom, start_pos, end_pos, 
                                                    truncate=True, max_depth=500):
                        depth = 0
                        for pileup_read in pileup_col.pileups:
                            if (not pileup_read.is_del and not pileup_read.is_refskip and
                                pileup_read.alignment.mapping_quality >= 8):
                                depth += 1
                        
                        region_depth += depth
                        region_length += 1
                    
                    if region_length > 0:
                        avg_depth = region_depth / region_length
                        if 0.5 <= avg_depth <= 500:  # Reasonable range
                            coverage_samples.append(avg_depth)
                    
                except Exception as e:
                    logger.debug(f"Error sampling {chrom}:{start_pos}-{end_pos}: {e}")
                    continue
        
        # Method 2: Quick read count estimation if pileup fails
        if len(coverage_samples) < 5:
            logger.info("Pileup method failed, trying read count method...")
            
            try:
                # Count reads in first chromosome
                if self.bam.nreferences > 0:
                    chrom = self.bam.references[0]
                    chrom_len = self.bam.lengths[0]
                    
                    # Sample 1Mb region
                    sample_len = min(1000000, chrom_len // 2)
                    sample_start = chrom_len // 4
                    sample_end = sample_start + sample_len
                    
                    read_count = 0
                    total_read_length = 0
                    
                    for read in self.bam.fetch(chrom, sample_start, sample_end):
                        if (not read.is_unmapped and not read.is_secondary and
                            not read.is_supplementary and read.mapping_quality >= 10):
                            read_count += 1
                            total_read_length += read.query_length or 100  # Default 100bp
                    
                    if read_count > 0:
                        # Estimate coverage: (reads * avg_read_length) / region_length
                        avg_read_length = total_read_length / read_count
                        estimated_coverage = (read_count * avg_read_length) / sample_len
                        
                        if 1 <= estimated_coverage <= 200:
                            coverage_samples = [estimated_coverage] * 10
                            logger.info(f"Estimated coverage from read count: {estimated_coverage:.1f}x")
            
            except Exception as e:
                logger.warning(f"Read count method failed: {e}")
        
        # Use samples to estimate coverage
        if len(coverage_samples) >= 5:
            # Remove outliers
            coverage_samples = sorted(coverage_samples)
            # Remove top and bottom 20%
            trim = len(coverage_samples) // 5
            if trim > 0:
                coverage_samples = coverage_samples[trim:-trim]
            
            median_coverage = np.median(coverage_samples)
            logger.info(f"Successfully estimated coverage: {median_coverage:.1f}x from {len(coverage_samples)} samples")
            return median_coverage
        
        # Final fallback: use expected coverage from params
        expected = self.expected_coverage
        logger.warning(f"Could not estimate coverage from BAM, using expected: {expected}x")
        logger.warning("This may affect genotyping and filtering accuracy")
        
        return expected 
        
    def run(self):
        logger.info("=== Complex_SV SV Detection Pipeline with Sniffles2-style Filtering ===")
        
        logger.info("1. Collecting all evidence...")
        self.collect_all_evidence()
        
        logger.info("2. Detecting deletions...")
        self.detect_deletions()
        
        logger.info("3. Detecting insertions...")
        self.detect_insertions()
        
        logger.info("4. Detecting ALL duplications...")
        self.detect_all_duplications_improved()
        
        logger.info("5. Detecting inversions...")
        self.detect_inversions()
        
        logger.info("6. Detecting BNDs...")
        self.detect_bnds_improved()
        
        logger.info("7. Detecting complex SVs (Sniffles2 method)...")
        self.detect_complex_svs_sniffles2_style()
        
        logger.info("8. Clustering evidence...")
        calls = self.cluster_evidence()
        
        logger.info("9. Identifying complex events from clustered calls...")
        calls = self.identify_complex_events(calls)
        
        logger.info("10. Genotyping...")
        calls = self.genotype_calls(calls)
        
        logger.info("11. Applying Sniffles2-style filters...")
        calls = self.apply_filters_sniffles2_style(calls)
        
        return calls
    
    def collect_all_evidence(self):
        """Collect all types of evidence"""
        total_reads = 0
        reads_with_evidence = 0
        
        for read in self.bam.fetch():
            total_reads += 1
            
            if read.is_unmapped or read.is_secondary:
                continue
            
            read_name = read.query_name
            has_evidence = False
            
            # Store segments for BND and complex SV detection
            segment = {
                'chr': read.reference_name,
                'start': read.reference_start,
                'end': read.reference_end,
                'strand': '-' if read.is_reverse else '+',
                'mapq': read.mapping_quality,
                'is_primary': not read.is_supplementary,
                'is_supplementary': read.is_supplementary,
                'cigar': read.cigarstring
            }
            self.segments[read_name].append(segment)
            
            # Process CIGAR for precise DEL/INS
            if not read.is_supplementary and read.mapping_quality >= self.min_mapq:
                if self.extract_cigar_evidence(read):
                    has_evidence = True
            
            # Process SA tag for all SV types
            if read.has_tag('SA'):
                if self.process_sa_tag(read):
                    has_evidence = True
            
            # Check for interspersed duplications
            if read.is_supplementary:
                self.check_for_interspersed_dup(read)
                has_evidence = True
            
            if has_evidence:
                reads_with_evidence += 1
        
        logger.info(f"  Processed {total_reads} reads, {reads_with_evidence} with SV evidence")
    
    def extract_cigar_evidence(self, read):
        """Extract DEL/INS from CIGAR and look for complex patterns"""
        ref_pos = read.reference_start
        cigar_events = []
        found_evidence = False
        
        for op, length in read.cigartuples:
            if op == 0 or op == 7 or op == 8:  # M/=/X
                ref_pos += length
            elif op == 1:  # Insertion
                if length >= self.min_sv_size:
                    evidence = Evidence('CIGAR_INS', read.reference_name, ref_pos, ref_pos,
                                      'INS', length, read.query_name, read.mapping_quality, {})
                    self.evidence['INS'].append(evidence)
                    cigar_events.append(('INS', ref_pos, length))
                    found_evidence = True
            elif op == 2:  # Deletion
                if length >= 15:
                    evidence = Evidence('CIGAR_DEL', read.reference_name, ref_pos, ref_pos + length,
                                      'DEL', length, read.query_name, read.mapping_quality, {})
                    self.evidence['DEL'].append(evidence)
                    cigar_events.append(('DEL', ref_pos, length))
                    found_evidence = True
                ref_pos += length
            elif op == 3:  # N
                ref_pos += length
        
        # Check for complex patterns
        if len(cigar_events) >= 2:
            self.check_for_complex_pattern(read, cigar_events)
        
        return found_evidence
    
    def check_for_complex_pattern(self, read, cigar_events):
        """Check if multiple CIGAR events form a complex SV"""
        if len(cigar_events) < 2:
            return
        
        for i in range(len(cigar_events) - 1):
            event1 = cigar_events[i]
            event2 = cigar_events[i + 1]
            
            distance = abs(event2[1] - event1[1])
            
            if distance <= self.complex_max_distance:
                extra = {
                    'components': [event1[0], event2[0]],
                    'positions': [event1[1], event2[1]],
                    'sizes': [event1[2], event2[2]]
                }
                
                evidence = Evidence('COMPLEX_CIGAR', read.reference_name, 
                                  min(event1[1], event2[1]),
                                  max(event1[1] + event1[2], event2[1] + event2[2]),
                                  'CPX', 0, read.query_name, read.mapping_quality, extra)
                self.complex_signatures.append(evidence)
    
    def process_sa_tag(self, read):
        """Process SA tag for split read evidence"""
        sa_tag = read.get_tag('SA')
        primary_chr = read.reference_name
        primary_start = read.reference_start
        primary_end = read.reference_end
        primary_strand = '-' if read.is_reverse else '+'
        
        sa_alignments = []
        found_evidence = False
        
        for sa_entry in sa_tag.rstrip(';').split(';'):
            if not sa_entry:
                continue
            
            parts = sa_entry.split(',')
            if len(parts) < 6:
                continue
            
            sa_chr = parts[0]
            sa_pos = int(parts[1]) - 1
            sa_strand = parts[2]
            sa_mapq = int(parts[4])
            
            if sa_mapq < self.min_mapq:
                continue
            
            sa_alignments.append({
                'chr': sa_chr,
                'pos': sa_pos,
                'strand': sa_strand,
                'mapq': sa_mapq
            })
            
            # BND - different chromosome or very distant
            if sa_chr != primary_chr:
                extra = {'chr2': sa_chr, 'pos2': sa_pos, 'strand2': sa_strand}
                # FIX: BND should have start=end (single breakpoint)
                evidence = Evidence('SA_BND', primary_chr, primary_end, primary_end,
                                  'BND', 0, read.query_name, min(read.mapping_quality, sa_mapq), extra)
                self.evidence['BND'].append(evidence)
                self.bnd_evidence.append(evidence)
                found_evidence = True
            
            # Inversion - different strand
            elif sa_strand != primary_strand:
                if abs(sa_pos - primary_end) >= self.min_sv_size:
                    evidence = Evidence('SA_INV', primary_chr, min(primary_start, sa_pos),
                                      max(primary_end, sa_pos), 'INV',
                                      abs(sa_pos - primary_start), read.query_name,
                                      min(read.mapping_quality, sa_mapq), {})
                    self.evidence['INV'].append(evidence)
                    found_evidence = True
            
            # Potential duplication - large gap on same strand
            elif abs(sa_pos - primary_end) > 5000:
                extra = {'target_pos': sa_pos, 'source_pos': primary_start}
                evidence = Evidence('SA_DUP', primary_chr, min(primary_start, sa_pos),
                                  max(primary_end, sa_pos), 'DUP',
                                  primary_end - primary_start, read.query_name,
                                  min(read.mapping_quality, sa_mapq), extra)
                self.evidence['DUP'].append(evidence)
                found_evidence = True
                
                if abs(sa_pos - primary_end) > 50000:
                    self.interspersed_evidence.append(evidence)
        
        # Check for complex patterns
        if len(sa_alignments) >= 2:
            self.check_for_complex_sa_pattern(read, primary_chr, primary_start, primary_end, sa_alignments)
        
        return found_evidence
    
    def check_for_complex_sa_pattern(self, read, primary_chr, primary_start, primary_end, sa_alignments):
        """Check if multiple SA alignments indicate complex SV"""
        components = []
        
        for sa in sa_alignments:
            if sa['chr'] != primary_chr:
                components.append('BND')
            elif sa['strand'] != ('+' if read.is_reverse else '-'):
                components.append('INV')
            else:
                distance = abs(sa['pos'] - primary_end)
                if distance > 50000:
                    components.append('DUP')
                elif distance > 1500:
                    components.append('DEL')
        
        if len(set(components)) >= 2:
            extra = {
                'components': components,
                'sa_count': len(sa_alignments)
            }
            
            evidence = Evidence('COMPLEX_SA', primary_chr, primary_start, primary_end,
                              'CPX', primary_end - primary_start, read.query_name,
                              read.mapping_quality, extra)
            self.complex_signatures.append(evidence)
    
    def check_for_interspersed_dup(self, read):
        """Check supplementary alignments for interspersed duplication signature"""
        primary_segs = [s for s in self.segments[read.query_name] if s['is_primary']]
        
        if not primary_segs:
            return
        
        primary = primary_segs[0]
        
        if primary['chr'] == read.reference_name:
            distance = abs(read.reference_start - primary['end'])
            
            if distance > 50000:
                extra = {
                    'source_start': primary['start'],
                    'source_end': primary['end'],
                    'target_start': read.reference_start,
                    'target_end': read.reference_end,
                    'distance': distance
                }
                evidence = Evidence('INTER_DUP', read.reference_name,
                                  min(primary['start'], read.reference_start),
                                  max(primary['end'], read.reference_end),
                                  'DUP', primary['end'] - primary['start'],
                                  read.query_name, read.mapping_quality, extra)
                self.evidence['DUP'].append(evidence)
                self.interspersed_evidence.append(evidence)
    
    def detect_complex_svs_sniffles2_style(self):
        """Detect complex SVs using Sniffles2's approach"""
        logger.info("  Detecting complex SVs (Sniffles2 method)...")
        
        # Group by chromosome and position
        chr_evidence = defaultdict(list)
        for sv_type, evidences in self.evidence.items():
            for ev in evidences:
                chr_evidence[ev.chrom].append((ev, sv_type))
        
        for chrom, evidence_list in chr_evidence.items():
            evidence_list.sort(key=lambda x: x[0].start)
            
            # Sniffles2 uses 1kb window for complex events
            window_size = 1000
            
            i = 0
            while i < len(evidence_list):
                window_evidence = [(evidence_list[i][0], evidence_list[i][1])]
                j = i + 1
                
                while j < len(evidence_list):
                    if evidence_list[j][0].start - evidence_list[i][0].start <= window_size:
                        window_evidence.append((evidence_list[j][0], evidence_list[j][1]))
                        j += 1
                    else:
                        break
                
                # Check if this forms a complex event (Sniffles2 criteria)
                if len(window_evidence) >= 2:
                    # Check for shared reads
                    read_sets = []
                    for ev, _ in window_evidence:
                        if hasattr(ev, 'read_name') and ev.read_name:
                            if ',' in ev.read_name:
                                read_sets.append(set(ev.read_name.split(',')))
                            else:
                                read_sets.append({ev.read_name})
                    
                    if read_sets:
                        shared_reads = set.intersection(*read_sets)
                        total_reads = set.union(*read_sets)
                        
                        if total_reads and len(shared_reads) / len(total_reads) > 0.3:
                            # This is a complex event
                            components = [ev[1] for ev in window_evidence]
                            
                            start = min(ev[0].start for ev in window_evidence)
                            end = max(ev[0].end for ev in window_evidence)
                            
                            evidence = Evidence(
                                'COMPLEX_SNIFFLES2',
                                chrom,
                                start,
                                end,
                                'CPX',
                                end - start,
                                ','.join(shared_reads),
                                int(np.mean([ev[0].mapq for ev in window_evidence])),
                                {'components': components, 'num_components': len(components)}
                            )
                            
                            self.evidence['CPX'].append(evidence)
                
                i = j if j > i + 1 else i + 1
        
        logger.info(f"    Found {len(self.evidence['CPX'])} complex SV candidates")
    
    def detect_deletions(self):
        logger.info(f"  Found {len(self.evidence['DEL'])} deletion signatures")
    
    def detect_insertions(self):
        logger.info(f"  Found {len(self.evidence['INS'])} insertion signatures")
    
    def detect_all_duplications_improved(self):
        """Improved duplication detection for both tandem and interspersed"""
        logger.info("  Running improved duplication detection...")
        
        # Depth-based detection
        window_size = 1500
        step_size = 900
        
        # Get DUP-specific thresholds from params
        dup_depth_threshold = self.params.get('dup_depth_threshold', 3.0)
        
        for chrom in self.bam.references[:min(5, len(self.bam.references))]:
            chrom_len = self.bam.get_reference_length(chrom)
            
            for start in range(0, min(chrom_len, 2000000), step_size):
                end = min(start + window_size, chrom_len)
                
                depth = 0
                reads_with_sa = []
                
                try:
                    for read in self.bam.fetch(chrom, start, end):
                        if not read.is_secondary and not read.is_supplementary:
                            if read.mapping_quality >= self.min_mapq:
                                depth += 1
                                
                                if read.has_tag('SA'):
                                    reads_with_sa.append(read)
                except:
                    continue
                
                expected_depth = self.coverage * window_size / 1000
                
                if expected_depth > 0 and depth > expected_depth * (dup_depth_threshold - 0.5):
                    distant_mappings = 0
                    for read in reads_with_sa:
                        sa_tag = read.get_tag('SA')
                        for sa in sa_tag.split(';'):
                            if sa and chrom in sa:
                                parts = sa.split(',')
                                if len(parts) > 1:
                                    try:
                                        sa_pos = int(parts[1])
                                        if abs(sa_pos - start) > 50000:
                                            distant_mappings += 1
                                    except:
                                        pass
                    
                    if distant_mappings >= 3:
                        dup_type = 'INTERSPERSED'
                        extra = {'distant_mappings': distant_mappings}
                    else:
                        dup_type = 'TANDEM'
                        extra = {}
                    
                    evidence = Evidence(f'DEPTH_{dup_type}', chrom, start, end,
                                      'DUP', end - start, f'depth_{chrom}_{start}', 30, extra)
                    self.evidence['DUP'].append(evidence)
                    
                    if dup_type == 'INTERSPERSED':
                        self.interspersed_evidence.append(evidence)
        
        tandem_count = len([e for e in self.evidence['DUP'] if 'TANDEM' in e.type])
        interspersed_count = len(self.interspersed_evidence)
        
        logger.info(f"  Found {len(self.evidence['DUP'])} total duplications")
        logger.info(f"    - Tandem: {tandem_count}")
        logger.info(f"    - Interspersed: {interspersed_count}")
    
    def detect_inversions(self):
        logger.info(f"  Found {len(self.evidence['INV'])} inversion signatures")
    
    def detect_bnds_improved(self):
        """Improved BND detection using segment analysis - FIXED POSITIONS"""
        logger.info("  Detecting BNDs with improved method...")
        
        for read_name, segments in self.segments.items():
            if len(segments) < 2:
                continue
            
            segments.sort(key=lambda x: x['start'])
            
            for i in range(len(segments) - 1):
                seg1 = segments[i]
                seg2 = segments[i + 1]
                
                if seg1.get('is_supplementary') and seg2.get('is_supplementary'):
                    continue
                
                # Inter-chromosomal BND
                if seg1['chr'] != seg2['chr']:
                    extra = {
                        'chr1': seg1['chr'],
                        'pos1': seg1['end'],
                        'chr2': seg2['chr'],
                        'pos2': seg2['start'],
                        'strand1': seg1['strand'],
                        'strand2': seg2['strand']
                    }
                    
                    if seg1['strand'] == '+' and seg2['strand'] == '+':
                        bnd_type = 'N['
                    elif seg1['strand'] == '+' and seg2['strand'] == '-':
                        bnd_type = 'N]'
                    elif seg1['strand'] == '-' and seg2['strand'] == '+':
                        bnd_type = '[N'
                    else:
                        bnd_type = ']N'
                    
                    extra['bnd_type'] = bnd_type
                    
                    # FIX: For BND, start=end (it's a breakpoint)
                    evidence = Evidence('SEGMENT_BND', seg1['chr'], seg1['end'], seg1['end'],
                                      'BND', 0, read_name, min(seg1['mapq'], seg2['mapq']), extra)
                    self.evidence['BND'].append(evidence)
                    self.bnd_evidence.append(evidence)
                
                # Intra-chromosomal distant BND
                elif abs(seg2['start'] - seg1['end']) > 100000:
                    extra = {
                        'chr1': seg1['chr'],
                        'pos1': seg1['end'],
                        'chr2': seg2['chr'],
                        'pos2': seg2['start'],
                        'distance': abs(seg2['start'] - seg1['end'])
                    }
                    # FIX: For BND, start=end
                    evidence = Evidence('DISTANT_BND', seg1['chr'], seg1['end'], seg1['end'],
                                      'BND', 0, read_name, min(seg1['mapq'], seg2['mapq']), extra)
                    self.evidence['BND'].append(evidence)
                    self.bnd_evidence.append(evidence)
        
        logger.info(f"  Found {len(self.bnd_evidence)} BND signatures")
    
    def cluster_evidence(self):
        """Cluster evidence with type-specific requirements"""
        calls = []
        
        for sv_type, evidences in self.evidence.items():
            if not evidences:
                continue
            
            chr_evidence = defaultdict(list)
            for ev in evidences:
                chr_evidence[ev.chrom].append(ev)
            
            for chrom, chrom_evs in chr_evidence.items():
                chrom_evs.sort(key=lambda x: x.start)
                
                clusters = []
                current_cluster = [chrom_evs[0]]
                
                for ev in chrom_evs[1:]:
                    # Type-specific clustering distances
                    if sv_type == 'BND':
                        max_dist = 900  # More permissive for BND
                    elif sv_type == 'DEL':
                        max_dist = 1000   # Moderate for DEL
                    elif sv_type in ['DUP']:
                        max_dist = 250
                    elif sv_type in ['INS']:
                        max_dist = 250   # Stricter for INS/DUP
                    elif sv_type == 'CPX':
                        max_dist = 70   # Strict for CPX
                    else:
                        max_dist = 200
                    
                    # Size similarity check
                    if sv_type not in ['BND', 'CPX']:
                        size_similar = abs(ev.size - current_cluster[-1].size) / max(ev.size, current_cluster[-1].size, 1) < 0.5
                    else:
                        size_similar = True
                    
                    if abs(ev.start - current_cluster[-1].start) <= max_dist and size_similar:
                        current_cluster.append(ev)
                    else:
                        if sv_type == 'BND':
                            min_support = 1  # Most permissive
                        elif sv_type == 'DEL':
                            min_support = 1  # Permissive for deletions
                        elif sv_type == 'DUP':
                            min_support = 3  # Moderate for duplications
                        elif sv_type == 'INS':
                            min_support = 2  # Strict for insertions
                        elif sv_type == 'CPX':
                            min_support = 7
                        
                        if len(current_cluster) >= min_support:
                            clusters.append(current_cluster)
                        current_cluster = [ev]
                
                # Handle last cluster with same type-specific requirements
                if sv_type == 'BND':
                    min_support = 1
                elif sv_type == 'DEL':
                    min_support = 1
                elif sv_type in ['DUP']:
                    min_support = 3
                elif sv_type in ['INS']:
                    min_support = 1
                elif sv_type == 'CPX':
                    min_support = 7
                else:
                    min_support = 2
                    
                if len(current_cluster) >= min_support:
                    clusters.append(current_cluster)
                
                for cluster in clusters:
                    call = self.create_call(cluster, sv_type)
                    if call:
                        calls.append(call)
        
        return calls
    
    def create_call(self, cluster, sv_type):
        """Create SV call from evidence cluster - FIXED END POSITION"""
        reads = set()
        depth_evidence_count = 0
        has_depth_evidence = False
        
        for ev in cluster:
            if ev.read_name.startswith('depth_'):
                depth_evidence_count += 5  # Balanced weight for depth evidence
                has_depth_evidence = True
            elif ev.read_name.startswith('complex_'):
                reads.add(ev.read_name)
            else:
                reads.add(ev.read_name)
        
        support = len(reads) + depth_evidence_count
        
        # Enhanced support for DUP
        if sv_type == 'DUP':
            sa_evidence = sum(1 for ev in cluster if 'SA_' in ev.type)
            depth_evidence = sum(1 for ev in cluster if 'DEPTH' in ev.type)
            
            if depth_evidence > 0:
                support = max(support, depth_evidence * 3.5 + sa_evidence)
            
            if has_depth_evidence and support < 4.5:
                support = 4.5
        
        # Enhanced support for BND
        if sv_type == 'BND':
            support = max(len(reads), len(cluster))
            if len(cluster) >= 2:
                support = max(support, len(cluster) * 2)
        
        # Handle complex SVs
        if sv_type == 'CPX':
            components = []
            for ev in cluster:
                if ev.extra and 'components' in ev.extra:
                    components.extend(ev.extra['components'])
            
            sig = SVSignature(
                chrom=cluster[0].chrom,
                start=int(np.median([ev.start for ev in cluster])),
                end=int(np.median([ev.end for ev in cluster])),
                sv_type='CPX',
                size=int(np.median([ev.size for ev in cluster])),
                support=support if support > 0 else len(cluster),
                mapq=int(np.mean([ev.mapq for ev in cluster])),
                reads=reads if reads else {f"complex_{cluster[0].chrom}_{cluster[0].start}"},
                af=support / max(1, self.coverage),
                is_complex=True,
                complex_components=components if components else ['Unknown'],
                source='complex_detection'
            )
            return sig
        
        # Type-specific minimum support
        if sv_type == 'BND':
            min_support = 1
        elif sv_type == 'DUP':
            has_depth = any('DEPTH' in ev.type for ev in cluster)
            min_support = 1 if has_depth else 2
        else:
            min_support = 2
        
        if support < min_support:
            return None
        
        starts = [ev.start for ev in cluster]
        ends = [ev.end for ev in cluster]
        
        # FIX: Ensure start < end for non-BND events
        median_start = int(np.median(starts))
        median_end = int(np.median(ends))
        
        # For BND, end should equal start (it's a breakpoint)
        if sv_type == 'BND':
            median_end = median_start  # BND is a single breakpoint
        else:
            # Ensure end > start for all other SV types
            if median_end <= median_start:
                # Calculate from size if available
                median_size = int(np.median([ev.size for ev in cluster]))
                if median_size > 0:
                    median_end = median_start + median_size
                else:
                    median_end = median_start + 1  # Minimum 1bp
        
        sig = SVSignature(
            chrom=cluster[0].chrom,
            start=median_start,
            end=median_end,
            sv_type=sv_type,
            size=int(np.median([ev.size for ev in cluster])) if sv_type != 'BND' else 0,
            support=support,
            mapq=int(np.mean([ev.mapq for ev in cluster])),
            reads=reads,
            af=support / max(1, self.coverage),
            std_pos=np.std(starts) if len(starts) > 1 else 0,
            std_len=np.std([ev.size for ev in cluster]) if len(cluster) > 1 and sv_type != 'BND' else 0,
            precise=any('CIGAR' in ev.type for ev in cluster),
            source='depth' if any('DEPTH' in ev.type for ev in cluster) else 'split_read'
        )
        
        # Check repeat context
        repeat_info = self.repeat_handler.check_repeat_context(sig.chrom, sig.start, sig.end)
        sig.in_repeat = repeat_info['in_repeat']
        sig.repeat_type = repeat_info['repeat_type']
        
        # Type-specific processing
        if sv_type == 'DUP':
            is_interspersed = any(ev in self.interspersed_evidence for ev in cluster)
            
            if is_interspersed or any('INTERSPERSED' in ev.type or 'INTER_' in ev.type for ev in cluster):
                sig.complex_type = 'DUP:INTERSPERSED'
            else:
                sig.complex_type = 'DUP:TANDEM'
            
            sig.cn = 2.0
            
            if any('DEPTH' in ev.type for ev in cluster):
                sig.source = 'depth'
                sig.depth_support = sum(1 for ev in cluster if 'DEPTH' in ev.type)
            if any('SA_' in ev.type for ev in cluster):
                sig.split_support = sum(1 for ev in cluster if 'SA_' in ev.type)
            
            # Boost AF for depth-based DUPs
            if sig.source == 'depth':
                sig.af = max(0.15, sig.af)
        
        elif sv_type == 'BND':
            for ev in cluster:
                if ev.extra and 'chr2' in ev.extra:
                    sig.mate_chr = ev.extra['chr2']
                    sig.mate_pos = ev.extra.get('pos2', sig.end)
                    sig.bnd_type = ev.extra.get('bnd_type', 'N[')
                    break
            
            if not sig.mate_chr:
                sig.mate_chr = sig.chrom
                sig.mate_pos = sig.end
                sig.bnd_type = 'N['
            
            sig.split_support = sum(1 for ev in cluster if 'SA_BND' in ev.type or 'SEGMENT_BND' in ev.type)
        
        elif sv_type == 'INV':
            sig.split_support = sum(1 for ev in cluster if 'SA_INV' in ev.type)
        
        elif sv_type in ['DEL', 'INS']:
            sig.split_support = sum(1 for ev in cluster if 'CIGAR' in ev.type)
            sig.spanning_support = len(reads) - sig.split_support
        
        return sig
    
    def identify_complex_events(self, calls):
        """Identify complex events from nearby simple SVs"""
        logger.info("  Identifying complex events from clustered calls...")
        
        chr_calls = defaultdict(list)
        for call in calls:
            chr_calls[call.chrom].append(call)
        
        complex_calls = []
        processed_indices = set()
        
        for chrom, chrom_calls in chr_calls.items():
            chrom_calls.sort(key=lambda x: x.start)
            
            for i in range(len(chrom_calls)):
                if i in processed_indices:
                    continue
                
                call1 = chrom_calls[i]
                nearby_calls = [call1]
                nearby_indices = {i}
                
                for j in range(i + 1, len(chrom_calls)):
                    call2 = chrom_calls[j]
                    
                    distance = call2.start - call1.end
                    
                    if distance <= self.complex_max_distance:
                        shared_reads = call1.reads.intersection(call2.reads)
                        
                        if shared_reads or distance <= self.complex_merge_distance:
                            nearby_calls.append(call2)
                            nearby_indices.add(j)
                    elif distance > self.complex_max_distance * 2:
                        break
                
                if len(nearby_calls) >= 2:
                    all_reads = set()
                    for call in nearby_calls:
                        all_reads.update(call.reads)
                    
                    read_sv_count = defaultdict(int)
                    for call in nearby_calls:
                        for read in call.reads:
                            read_sv_count[read] += 1
                    
                    multi_sv_reads = sum(1 for count in read_sv_count.values() if count >= 2)
                    
                    if multi_sv_reads >= self.min_support_complex:
                        complex_call = SVSignature(
                            chrom=chrom,
                            start=min(call.start for call in nearby_calls),
                            end=max(call.end for call in nearby_calls),
                            sv_type='CPX',
                            size=max(call.end for call in nearby_calls) - min(call.start for call in nearby_calls),
                            support=len(all_reads),
                            mapq=int(np.mean([call.mapq for call in nearby_calls])),
                            reads=all_reads,
                            is_complex=True,
                            complex_components=[call.sv_type for call in nearby_calls],
                            af=len(all_reads) / max(1, self.coverage)
                        )
                        
                        complex_calls.append(complex_call)
                        processed_indices.update(nearby_indices)
                        
                        logger.info(f"    Found complex SV at {chrom}:{complex_call.start}-{complex_call.end} "
                                  f"with components: {complex_call.complex_components}")
        
        if complex_calls:
            logger.info(f"    Identified {len(complex_calls)} complex SVs from simple SV clustering")
            calls.extend(complex_calls)
        
        return calls
    
    def genotype_calls(self, calls):
        """Genotype SV calls - WITH POSITION VALIDATION"""
        for call in calls:
            # FIX: Validate and fix positions
            if call.sv_type == 'BND':
                # BND should have start = end
                call.end = call.start
                fetch_start = max(0, call.start - 1000)
                fetch_end = call.start + 1000
            else:
                # Ensure start < end for non-BND
                if call.end <= call.start:
                    if call.size > 0:
                        call.end = call.start + call.size
                    else:
                        call.end = call.start + 1
                
                fetch_start = max(0, call.start - 1000)
                fetch_end = call.end + 1000
            
            # Ensure fetch coordinates are valid
            if call.chrom in self.ref_lengths:
                chrom_len = self.ref_lengths[call.chrom]
                fetch_end = min(fetch_end, chrom_len)
                fetch_start = min(fetch_start, chrom_len - 1)
            
            # Get reads in region
            region_reads = set()
            try:
                for read in self.bam.fetch(call.chrom, fetch_start, fetch_end):
                    if not read.is_secondary and not read.is_supplementary:
                        if read.mapping_quality >= self.min_mapq:
                            region_reads.add(read.query_name)
            except ValueError as e:
                logger.warning(f"Could not fetch reads for {call.chrom}:{fetch_start}-{fetch_end}: {e}")
                call.genotype = './.'
                call.genotype_qual = 0
                call.alt_support = 0
                call.ref_support = 0
                continue
            
            sv_reads = call.reads
            ref_reads = region_reads - sv_reads
            
            call.alt_support = len(sv_reads)
            call.ref_support = len(ref_reads)
            total = call.alt_support + call.ref_support
            
            if total == 0:
                call.genotype = './.'
                call.genotype_qual = 0
            else:
                af = call.alt_support / total
                call.af = af
                
                # IMPROVED: More lenient genotyping thresholds
                if af < 0.08:  # Reduced from 0.15
                    call.genotype = '0/0'
                elif af < 0.55:  # Reduced from 0.65
                    call.genotype = '0/1'
                else:
                    call.genotype = '1/1'
                
                # BOOST AF for depth-based calls
                if hasattr(call, 'source') and call.source == 'depth':
                    call.af = max(call.af, 0.25)  # Minimum AF for depth calls
                
                call.genotype_qual = min(99, int(40 * af))  # Increased from 30
            
        return calls
    
    def calculate_strand_support(self, call):
        """Calculate strand support for a call"""
        strand_support = {'forward': 0, 'reverse': 0}
        
        try:
            for read_name in list(call.reads)[:50]:  # Sample first 100 reads
                for read in self.bam.fetch(call.chrom, max(0, call.start - 1000), call.end + 1000):
                    if read.query_name == read_name:
                        if read.is_reverse:
                            strand_support['reverse'] += 1
                        else:
                            strand_support['forward'] += 1
                        break
        except:
            # If can't calculate, assume balanced
            half = call.support // 2
            strand_support = {'forward': half, 'reverse': call.support - half}
        
        return strand_support
            
    def apply_filters_sniffles2_style(self, calls):
        """Apply Sniffles2-inspired filtering approach - BALANCED FOR DUPS"""
        
        for call in calls:
            if not call.strand_support:
                call.strand_support = self.calculate_strand_support(call)
            
            call.qual = self.quality_scorer.calculate_sv_quality(call)
            
            filters = []
            
            # Get base thresholds
            min_support_base = self.params.get(f'min_support_{call.sv_type.lower()}', 3)
            min_af_base = self.params.get(f'min_af_{call.sv_type.lower()}', 0.15)
            
            # Repeat-aware adjustments
            if call.in_repeat:
                if call.sv_type == 'CPX':
                    modifier = self.params.get('repeat_support_modifier_cpx', 0.7)
                    min_support = max(1, int(min_support_base * modifier))
                    min_af = min_af_base * 0.8
                elif call.sv_type in ['DEL', 'INS', 'INV']:
                    modifier = self.params.get(f'repeat_support_modifier_{call.sv_type.lower()}', 1.3)
                    min_support = int(min_support_base * modifier)
                    min_af = min_af_base * 1.1
                elif call.sv_type == 'DUP':
                    # BALANCED FOR DUPS in repeats
                    if call.source == 'depth' or call.depth_support > 2:
                        min_support = 2  # more lenient with depth support
                        min_af = self.params.get('dup_min_af_depth', 0.30)
                    else:
                        min_support = 3
                        min_af = self.params.get('dup_min_af_split', 0.30)
                else:
                    min_support = min_support_base
                    min_af = min_af_base
            else:
                # Non-repeat adjustments
                if call.sv_type == 'BND':
                    min_support = max(1, min_support_base)
                    min_af = min_af_base * 0.8
                elif call.sv_type == 'DEL':
                    min_support = max(1, min_support_base - 1) 
                    min_af = min_af_base * 0.4
                elif call.sv_type in ['INS', 'CPX']:
                    min_support = min_support_base + 1
                    min_af = min_af_base * 1.1
                elif call.sv_type == 'DUP':
                    # BALANCED FOR DUPS outside repeats
                    if call.source == 'depth' and call.depth_support >= 2:
                        min_support = 2
                        min_af = self.params.get('dup_min_af_depth', 0.40)
                    else:
                        min_support = 3
                        min_af = self.params.get('dup_min_af_split', 0.40)
                else:
                    min_support = min_support_base
                    min_af = min_af_base
            
            # Apply filters
            if call.support < min_support:
                filters.append('MinSupport')
            
            if call.af < min_af:
                filters.append('MinAF')
            
            # Quality filter (lenient for DUPs)
            if call.sv_type == 'DUP':
                if call.qual < 10:
                    filters.append('LowQual')
            else:
                qual_threshold = 4 if call.sv_type in ['BND', 'DEL'] else 10
                if call.qual < qual_threshold:
                    filters.append('LowQual')
            
            # Special rescue rules for DUPs (BALANCED)
            if call.sv_type == 'DUP' and len(filters) > 0:
                if call.source == 'depth' and call.depth_support >= 3 and call.qual >= 20:
                    # Keep only critical filters
                    filters_to_keep = []
                    if 'MinSize' in filters and call.size < 50:
                        filters_to_keep.append('MinSize')
                    if 'LowQual' in filters and call.qual < 5:
                        filters_to_keep.append('LowQual')
                    filters = filters_to_keep
                elif call.qual >= 20:
                    # Rescue high-quality DUPs
                    filters = [f for f in filters if f not in ['MinAF', 'ImpreciseBreakpoint']]
            
            # Rescue BNDs with any evidence
            if call.sv_type == 'BND' and call.support >= 1:
                filters = [f for f in filters if f != 'MinSupport']
            
            # Rescue CPX in repeats
            if call.sv_type == 'CPX' and call.in_repeat and call.support >= 2:
                filters = [f for f in filters if f not in ['MinSupport', 'MinAF']]
            
            call.filter = 'PASS' if not filters else ';'.join(filters)
        
        return calls


    def log_filter_statistics(self, calls):
        """Log filtering statistics (Sniffles2-style)"""
        total = len(calls)
        passed = sum(1 for c in calls if c.filter == 'PASS')
        
        logger.info("="*60)
        logger.info("FILTERING STATISTICS (Sniffles2-style):")
        logger.info(f"Total calls: {total}")
        logger.info(f"PASS calls: {passed} ({100*passed/max(1,total):.1f}%)")
        
        # Per-type statistics
        for sv_type in ['DEL', 'INS', 'DUP', 'INV', 'BND', 'CPX']:
            type_calls = [c for c in calls if c.sv_type == sv_type]
            type_pass = [c for c in type_calls if c.filter == 'PASS']
            
            if type_calls:
                avg_qual = np.mean([c.qual for c in type_calls])
                logger.info(f"  {sv_type}: {len(type_pass)}/{len(type_calls)} PASS "
                            f"({100*len(type_pass)/len(type_calls):.1f}%), "
                            f"Avg QUAL: {avg_qual:.1f}")
        
        # Filter reason statistics
        filter_counts = defaultdict(int)
        for call in calls:
            if call.filter != 'PASS':
                for f in call.filter.split(';'):
                    filter_counts[f] += 1
        
        if filter_counts:
            logger.info("\nFilter reasons:")
            for filter_name, count in sorted(filter_counts.items(), key=lambda x: x[1], reverse=True):
                logger.info(f"  {filter_name}: {count} ({100*count/total:.1f}%)")
        
        # Repeat region statistics
        in_repeat = sum(1 for c in calls if c.in_repeat)
        if in_repeat > 0:
            logger.info(f"\nCalls in repeat regions: {in_repeat} ({100*in_repeat/total:.1f}%)")
            repeat_pass = sum(1 for c in calls if c.in_repeat and c.filter == 'PASS')
            logger.info(f"  PASS in repeats: {repeat_pass} ({100*repeat_pass/max(1,in_repeat):.1f}%)")
        
        logger.info("="*60)

    def write_vcf(self, calls, output_file):
        """Write VCF with all proper headers - FIXED END POSITION"""
        with open(output_file, 'w') as f:
            # Header
            f.write("##fileformat=VCFv4.3\n")
            f.write("##source=Complex_SVe\n")
            
            # Contig headers
            for chrom, length in self.ref_lengths.items():
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
            f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
            f.write('##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Position StDev">\n')
            f.write('##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Length StDev">\n')
            f.write('##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise SV">\n')
            f.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise SV">\n')
            f.write('##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Mean MAPQ">\n')
            f.write('##INFO=<ID=DUPTYPE,Number=1,Type=String,Description="Duplication type">\n')
            f.write('##INFO=<ID=CN,Number=1,Type=Float,Description="Copy number">\n')
            f.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromosome">\n')
            f.write('##INFO=<ID=POS2,Number=1,Type=Integer,Description="Mate position">\n')
            f.write('##INFO=<ID=CPX_TYPE,Number=.,Type=String,Description="Complex SV components">\n')
            f.write('##INFO=<ID=REPEAT,Number=1,Type=String,Description="Repeat type if in repeat region">\n')
            f.write('##INFO=<ID=STRAND_RATIO,Number=1,Type=Float,Description="Strand bias ratio">\n')
            
            # FILTER fields
            f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
            f.write('##FILTER=<ID=MinSize,Description="SV size below minimum threshold">\n')
            f.write('##FILTER=<ID=MinSupport,Description="Read support below minimum threshold">\n')
            f.write('##FILTER=<ID=MinAF,Description="Allele frequency below minimum threshold">\n')
            f.write('##FILTER=<ID=LowMAPQ,Description="Low mapping quality">\n')
            f.write('##FILTER=<ID=LowQual,Description="Low quality score">\n')
            f.write('##FILTER=<ID=LowGQ,Description="Low genotype quality">\n')
            f.write('##FILTER=<ID=HighSTDPos,Description="High standard deviation of breakpoint position">\n')
            f.write('##FILTER=<ID=StrandBias,Description="Extreme strand bias">\n')
            f.write('##FILTER=<ID=Artifact,Description="Likely artifact">\n')
            f.write('##FILTER=<ID=ImpreciseBreakpoint,Description="Imprecise breakpoint">\n')
            
            # FORMAT fields
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
            f.write('##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Ref support">\n')
            f.write('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Alt support">\n')
            
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        
            # Write variants
            for i, call in enumerate(sorted(calls, key=lambda x: (x.chrom, x.start))):
                chrom = call.chrom
                pos = call.start + 1  # VCF is 1-based
                sv_id = f"{call.sv_type}_{i+1}"
                
                # FIX: Ensure END is valid
                if call.sv_type == 'BND':
                    # For BND, END should be same as POS
                    end_pos = pos
                else:
                    # For other SVs, END must be > POS
                    end_pos = call.end + 1
                    if end_pos <= pos:
                        # Fix invalid END
                        if call.size > 0:
                            end_pos = pos + call.size
                        else:
                            end_pos = pos + 1
                
                # REF and ALT
                ref = 'N'
                if call.sv_type == 'CPX':
                    alt = '<CPX>'
                elif call.sv_type == 'BND':
                    if call.mate_chr and call.mate_pos is not None:
                        # FIX: Ensure mate_pos is valid
                        mate_pos = max(1, call.mate_pos + 1)
                        if call.bnd_type == 'N[':
                            alt = f"N[{call.mate_chr}:{mate_pos}["
                        elif call.bnd_type == 'N]':
                            alt = f"N]{call.mate_chr}:{mate_pos}]"
                        elif call.bnd_type == '[N':
                            alt = f"[{call.mate_chr}:{mate_pos}[N"
                        else:
                            alt = f"]{call.mate_chr}:{mate_pos}]N"
                    else:
                        alt = '<BND>'
                elif call.sv_type == 'DUP' and call.complex_type:
                    alt = f"<{call.complex_type}>"
                else:
                    alt = f"<{call.sv_type}>"
                
                qual = int(call.qual) if call.qual > 0 else 10
                
                # INFO - with fixed END
                info_parts = [
                    f"SVTYPE={call.sv_type}",
                    f"END={end_pos}",  # Use fixed end_pos
                    f"SVLEN={call.size if call.sv_type != 'DEL' else -call.size}",
                    f"SUPPORT={call.support}",
                    f"AF={call.af:.3f}",
                    f"STDEV_POS={call.std_pos:.1f}",
                    f"STDEV_LEN={call.std_len:.1f}",
                    f"MAPQ={call.mapq}"
                ]
                
                if call.precise:
                    info_parts.append("PRECISE")
                else:
                    info_parts.append("IMPRECISE")
                
                if call.sv_type == 'DUP':
                    dup_type = 'INTERSPERSED' if 'INTERSPERSED' in str(call.complex_type) else 'TANDEM'
                    info_parts.append(f"DUPTYPE={dup_type}")
                    info_parts.append(f"CN={call.cn:.1f}")
                
                if call.sv_type == 'BND' and call.mate_chr:
                    info_parts.append(f"CHR2={call.mate_chr}")
                    if call.mate_pos is not None:
                        # FIX: Ensure valid mate position
                        info_parts.append(f"POS2={max(1, call.mate_pos+1)}")
                
                if call.sv_type == 'CPX' and call.complex_components:
                    components_str = ','.join(call.complex_components)
                    info_parts.append(f"CPX_TYPE={components_str}")
                
                if call.in_repeat and call.repeat_type:
                    info_parts.append(f"REPEAT={call.repeat_type}")
                
                if call.strand_support:
                    forward = call.strand_support.get('forward', 0)
                    reverse = call.strand_support.get('reverse', 0)
                    total = forward + reverse
                    if total > 0:
                        ratio = min(forward, reverse) / total
                        info_parts.append(f"STRAND_RATIO={ratio:.3f}")
                
                info = ';'.join(info_parts)
                
                # FORMAT
                format_str = "GT:GQ:DR:DV"
                sample_str = f"{call.genotype}:{call.genotype_qual}:{call.ref_support}:{call.alt_support}"
                
                f.write(f"{chrom}\t{pos}\t{sv_id}\t{ref}\t{alt}\t{qual}\t{call.filter}\t{info}\t{format_str}\t{sample_str}\n")
        
        logger.info(f"Wrote {len(calls)} variants to {output_file}")

    def close(self):
        self.bam.close()
        if self.ref:
            self.ref.close()

def main():
    if len(sys.argv) < 5:
        print("Usage: enhanced_sv_caller.py <bam> <ref> <params.json> <output.vcf> [trf.bed]")
        sys.exit(1)

    trf_bed = sys.argv[5] if len(sys.argv) > 5 else None

    caller = EnhancedSVCaller(sys.argv[1], sys.argv[2], sys.argv[3], trf_bed)
    calls = caller.run()
    caller.write_vcf(calls, sys.argv[4])

    # Summary
    sv_counts = Counter(c.sv_type for c in calls)
    pass_calls = [c for c in calls if c.filter == 'PASS']

    logger.info("="*60)
    logger.info("FINAL RESULTS:")

    for sv_type in ['DEL', 'INS', 'DUP', 'INV', 'BND', 'CPX']:
        total = sv_counts.get(sv_type, 0)
        passed = len([c for c in pass_calls if c.sv_type == sv_type])
        logger.info(f"  {sv_type}: {total} total, {passed} PASS")
        
        if sv_type == 'DUP' and total > 0:
            tandem = len([c for c in calls if c.sv_type == 'DUP' and 'TANDEM' in str(c.complex_type)])
            interspersed = len([c for c in calls if c.sv_type == 'DUP' and 'INTERSPERSED' in str(c.complex_type)])
            logger.info(f"    - Tandem: {tandem}, Interspersed: {interspersed}")
        elif sv_type == 'CPX' and total > 0:
            component_counts = defaultdict(int)
            for call in calls:
                if call.sv_type == 'CPX' and call.complex_components:
                    for comp in call.complex_components:
                        component_counts[comp] += 1
            if component_counts:
                logger.info(f"    - Components: {dict(component_counts)}")

    logger.info(f"\nTotal: {len(calls)} calls, {len(pass_calls)} PASS")
    logger.info(f"Coverage: {caller.coverage:.1f}x")
    logger.info("="*60)            

    caller.close()

if __name__ == "__main__":
    main()
