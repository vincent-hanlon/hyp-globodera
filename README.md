This repository contains scripts relating to the article

**Diversity at the HYP1 locus in potato cyst nematodes does not result from developmentally-programmed somatic mutations**

by Vincent C. T. Hanlon, Unnati Sonawala, Cian A. A. Raza, Luisa Kalkert, George Harpum, Lukas A. Burkhardt, Johannes Helder, and Sebastian Eves-van den Akker

### HYP1 allele parsing analysis

### Alignment plots for manual examination of FASTQ reads

### Nucleotide diversity calculation

This simple script calculates the nucleotide diversity across genomic regions (BED file) from biallelic or multiallelic SNVs in an all-sites VCF (should include invariant sites). The VCF must contain a single sample originating from a pool of individuals, under the assumption that there are many more individuals than reads. This way, every read is likely to represent a distinct individual. The script uses the allele depth field as a proxy for allele frequencies.
