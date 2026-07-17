#!/usr/bin/env python3


# This script calculates the nucleotide diversity across genomic regions (BED file) from biallelic or multiallelic SNVs in an all-sites VCF (should include invariant sites).
# The VCF must contain a single sample originating from a pool of individuals, under the assumption that there are many more individuals than reads. This way,
# every read is likely to represent a distinct individual. The script uses the allele depth field as a proxy for allele frequencies.
#
# Usage is python nucleotide_diversity_pi_from_pooled_vcf.py -b file1.bed -v file2.vcf -o output_file.txt



import pysam
import sys
import os
import argparse

def validate_vcf(vcf_path):
    """Ensures the VCF is bgzipped and indexed."""
    if not (vcf_path.endswith('.gz') or vcf_path.endswith('.bgz')):
        raise ValueError(f"Input VCF '{vcf_path}' must be bgzipped (e.g., .vcf.gz)")
    
    if not (os.path.exists(vcf_path + '.tbi') or os.path.exists(vcf_path + '.csi')):
        raise FileNotFoundError(f"Index file not found for '{vcf_path}'. Please index with 'tabix -p vcf <file>'")


def main():
    parser = argparse.ArgumentParser(description="Calculate multiallelic pi for a pooled sample from an all-sites VCF.")
    parser.add_argument("-v", "--vcf", required=True, help="bgzipped and tabix-indexed all-sites VCF")
    parser.add_argument("-b", "--bed", required=True, help="BED file of genomic intervals")
    parser.add_argument("-o", "--out", required=True, help="Output BED file")
    args = parser.parse_args()

    validate_vcf(args.vcf)

    vcf_in = pysam.VariantFile(args.vcf)
    
    # Assume there is only one sample in the VCF (the pool)
    sample_name = list(vcf_in.header.samples)[0]

    with open(args.bed, 'r') as bed_in, open(args.out, 'w') as bed_out:
        for line in bed_in:
            if line.startswith("#") or line.strip() == "":
                continue
            
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])

            L = 0              # Total called sites
            pi_sum = 0.0       # Sum of per-site nucleotide diversity

            try:
                for rec in vcf_in.fetch(chrom, start, end):
                    if not (start <= rec.pos - 1 < end):
                        continue
                    
                    sample_data = rec.samples[sample_name]
                    ad = sample_data.get('AD')

                    # Skip sites with missing data
                    if ad is None:
                        L += 1
                        continue
                    
                    # Calculate effective depth safely strictly from AD to avoid GATK DP inflation
                    eff_dp = sum(count for count in ad if count is not None)
                    if eff_dp <= 1:
                        continue
                    
                    L += 1

                    # Multiallelic Nucleotide Diversity (Pi)
                    sum_p_squared = 0.0
                    
                    for count in ad:
                        if count is not None:
                            # Add to Pi calculation
                            p_i = count / eff_dp
                            sum_p_squared += p_i**2
                            
                    finite_correction = eff_dp / (eff_dp - 1)
                    pi_site = finite_correction * (1.0 - sum_p_squared)
                    pi_sum += pi_site
                    
            except ValueError:
                # In case chromosome from BED is missing in the VCF entirely
                pass

            # Calculate interval-wide statistics
            if L > 0:
                pi_window = pi_sum / L

            else:
                pi_window = 0.0

            # Output: chrom, start, end, Pi
            out_cols = cols[:3] + [f"{pi_window:.6e}"]
            bed_out.write('\t'.join(out_cols) + '\n')

if __name__ == "__main__":
    main()
