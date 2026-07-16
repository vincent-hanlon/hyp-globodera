This repository contains scripts relating to the article

**Diversity at the HYP1 locus in potato cyst nematodes does not result from developmentally-programmed somatic mutations**

by Vincent C. T. Hanlon, Unnati Sonawala, Cian A. A. Raza, Luisa Kalkert, George Harpum, Lukas A. Burkhardt, Johannes Helder, and Sebastian Eves-van den Akker

### HYP1 allele parsing analysis

Given a FASTQ file of reads, the script [parse_hyp_motifs.py](parse_hyp_motifs.py) identifies reads that contain HYPs by checking if they have a perfect match either for a known motif or for one of six 18mers adjacent to the repetitive domain in exon 2. Then, it attempts to locate and extract the repetitive domain (in a consistent orientation) by aligning two conserved sequences that flank it. Finally, it parses the repetitive domain by substituting in the known motif sequences, by default allowing "fuzzy matches" where 1 error in the motif sequence is permitted. It outputs a table containing the extracted repetitive domain for each read, a parsed version of the domain both with amino acids and with motif names, a string of PHRED base qualities, etc. 

It is run like this:

`python parse_hyp_motifs.py fastq_file conserved_primers_fasta motif_table`

Examples of motif tables and the conserved primer sequences can be found in the directory [parsing_inputs/](parsing_inputs)

Subsequently, the output file can be processed to identify putative rare alleles in several ways. We used the routine found in the script [parsed_post_processing.sh](parsed_post_processing.sh).

### Alignment plots for manual examination of FASTQ reads

### Nucleotide diversity calculation

This simple script calculates the nucleotide diversity across genomic regions (BED file) from biallelic or multiallelic SNVs in an all-sites VCF (should include invariant sites). The VCF must contain a single sample originating from a pool of individuals, under the assumption that there are many more individuals than reads. This way, every read is likely to represent a distinct individual. The script uses the allele depth field as a proxy for allele frequencies.
