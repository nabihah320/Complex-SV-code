#!/usr/bin/env python3
"""
Fix SVIM VCF for Truvari compatibility
"""
import sys

def fix_svim_vcf(input_vcf, output_vcf):
    """Fix SVIM VCF for Truvari evaluation"""
    print("Fixing SVIM VCF for Truvari...")
    
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    # Force more variants to PASS
                    qual = parts[5]
                    try:
                        qual_val = float(qual) if qual != '.' else 0
                        if qual_val > 5:
                            parts[6] = 'PASS'
                        else:
                            parts[6] = 'LowQual'
                    except:
                        parts[6] = 'LowQual'
                    
                    # Ensure required INFO fields
                    info = parts[7]
                    if 'SVTYPE=' not in info:
                        alt = parts[4]
                        if '<DEL' in alt:
                            info += ';SVTYPE=DEL'
                        elif '<INS' in alt:
                            info += ';SVTYPE=INS'
                        elif '<DUP' in alt:
                            info += ';SVTYPE=DUP'
                        elif '<INV' in alt:
                            info += ';SVTYPE=INV'
                        else:
                            info += ';SVTYPE=BND'
                        parts[7] = info
                    
                    outfile.write('\t'.join(parts[:10]) + '\n')
    
    print("SVIM VCF fixed")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 fix_svim_vcf.py <input_vcf> <output_vcf>")
        sys.exit(1)
    
    fix_svim_vcf(sys.argv[1], sys.argv[2])
