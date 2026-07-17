# This messy-looking command really just selects reads that could be completely substituted with motifs (no extra bases left over), but with a little bit of flexibility at the start or end of the repeitive domain to account for any alignment problems.

# selecting the parsed representations of all reads
cut -f3 parsed_fastq_input.txt | \
# removing reads that couldn't be parsed into motifs at all
 grep -v -E '^[AGCT]+$' | \
# ignoring up to five consecutive nucleotides not substituted by a motif at the start of the repetitive domain 
 sed -E 's/^[[:space:]]{0,5}[ACGTN]{0,5}[[:space:]]{1}//' | \
# the same, but for the end of the repetitive domain
 sed -E 's/[[:space:]]{1}[ACGTN]{0,5}[[:space:]]{0,5}$//'| \
# removing any reads that couldn't be completely substituted (i.e., after subsituting with motifs, some DNA was leftover in the repetitive domain). Also removing the header.
 grep -v -E ' [AGCT]+ ' | grep -v 'protein' | \
# removing reads with 6 or more un-substituted nucleotides at the start or end of the repetitive domain (see the sed lines above)
 grep -v -E '^[AGCT]{6,} ' | grep -v -E ' [AGCT]{6,}$' | \
# formatting...
 sed -E 's/[[:space:]]+/ /g' | \
# counting the number of occurrences of each parsed repetitive domain allele
 sort | uniq -c  | sort -nr  > parsed_fastq_input_unique_counts.txt
