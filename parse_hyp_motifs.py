# This script takes a fastq file of reads, selects those that contain HYP1, locates the repetitive domain, and attempts to parse
# the repetitive domain by substituting in motif names/sequences using a matching routine that allows 1 error per motif. The output
# file also contains a reoriented DNA sequence (if it was reverse-complemented) along with a PHRED quality string. This file requires 
# further analysis (see github repo) to identify unique rare alelles. 
#
# The user must provide a fastq of reads, a fasta with the sequences that flank the HYP repetitive domain, and a table of motifs with
# three columns: DNA sequence, amino acid sequence, motif name
# 
# The script is not particularly fast, so if it is run on a large dataset it may be worth splitting up the FASTQ file and running 
# in parallel on the pieces.
#
# Usage: python parse_hyp_motifs.py fastq_file conserved_primers_fasta motif_table

import sys

fasta_filename=str(sys.argv[2])
fastq_filename=str(sys.argv[1])
motif_filename=str(sys.argv[3])

from Bio import SeqIO, Align
import pandas as pd
import re
import os
from math import log
from math import *
import pyfastx as px
import csv
import regex

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_sequence = ''.join([complement_dict[base] for base in reversed(dna_sequence)])
    return reverse_complement_sequence

def find_hvd(before, after, read, phred, match_thresh=0.4):
    #This locates the hvd (now called the "repetitive domain") by aligning conserved sequences that are adjacent to it,  provided by the user
    aligner = Align.PairwiseAligner(match_score=1.0, mode="global", mismatch_score=-1.0, open_gap_score=-3.0, extend_gap_score=-0.5, end_gap_score=0)

    precede=aligner.align(before, read, strand='+')[0]
    precede_end=precede.coordinates[1][-2].tolist()
    follow=aligner.align(after,read, strand='+')[0]
    follow_start=follow.coordinates[1][1].tolist()
    
    if not (precede_end < follow_start and precede.score/len(before) >= match_thresh and follow.score/len(after) >= match_thresh):
        read=reverse_complement(read)
        phred=phred[::-1]
        precede=aligner.align(before,read, strand='+')[0]
        precede_end=precede.coordinates[1][-2].tolist()
        follow=aligner.align(after,read, strand='+')[0]
        follow_start=follow.coordinates[1][1].tolist()
        
        if not (precede_end < follow_start and precede.score/len(before) > match_thresh and follow.score/len(after) > match_thresh):
            return None

    hvd=[read[precede_end-20:follow_start + 20], phred[precede_end-20:follow_start+20]]
    
    return(hvd)

def fuzzy_replace_short_strings(long_string, mapping_file_path, which="protein", errors=1):
    #This does a basic find a replace, except that mismatches are allowed, to parse a DNA sequence into amino acid motifs
    # Read the tab-delimited file and process replacements
    with open(mapping_file_path, "r", encoding="utf-8") as f:
        for line in f:
            short_str, protein, name = line.strip().split("\t")
            
            # Case-insensitive match for any occurrence, including substrings
            pattern = re.escape(short_str)
            
            # First doing the replacement with perfect matches (those are best)
            # Replace with name flanked by tab spaces
            if which=="protein":
                replacement = f" {protein} "
            elif which=="name":
                replacement = f" {name} "
            else:
                raise ValueError("permitted values for which are 'name' and 'protein'")
         
            long_string = re.sub(pattern, replacement, long_string, flags=re.IGNORECASE)

        # Now we can accept error-tolerant substitutions (including short gaps)
        if errors>0:
            f.seek(0)
            for line in f:
                short_str, protein, name = line.strip().split("\t")

            # Case-insensitive match for any occurrence, including substrings
                pattern = rf"(?bi)({regex.escape(short_str)}){{e<=1}}"
            # Replace with name flanked by tab spaces
                if which=="protein":
                    replacement = f" {protein} "
                elif which=="name":
                    replacement = f" {name} "
                else:
                    raise ValueError("permitted values for which are 'name' and 'protein'")

                long_string=regex.sub(pattern, replacement, long_string)

    return long_string


def Phred_ASCII_to_prob(symbol, offset=33):
    prob=10**-((ord(symbol)-offset)/10.0)
    return(prob)


def mean_read_qual(QUAL, offset=33):
    # Given an ASCII string representing per-base qualities, this returns the mean Phred quality of the read (or partial read)
    probabilities=[Phred_ASCII_to_prob(base,offset) for base in QUAL]
    avg=sum(probabilities)/len(probabilities)
    phred=-10*log(avg,10)
    return(phred)



# reading in the conserved sequences that will help us locate the hvd ("repetitive domain")
fasta = [seq for name, seq in px.Fasta(fasta_filename, build_index=False)]
before=fasta[0]
after=fasta[1]

# breaking up the conserved sequences into 18mers
n=18
before_pieces = [before[len(before)-i*n-n:len(before)-i*n] for i in range(0, 3)]
after_pieces = [after[i*n:i*n + n] for i in range(0, 3)]

# reading in the sequence motifs andf their revere complements 
patterns=pd.read_csv(motif_filename, sep="\t", header=None).iloc[:,0].tolist() + before_pieces + after_pieces #+ [before[len(before)-18:len(before)], after[0:18]]
rev_patterns=[reverse_complement(p) for p in patterns]
patterns=patterns + rev_patterns 

# reading a FASTQ of reads that might have HYPs
parsedlist=[]
fq = px.Fastq(fastq_filename, build_index=False)

for name, seq, qual in fq:
    # selecting reads that have at least one exact match of either a motif or of an 18mer from the nearby conserved sequences---these are probably HYP1s 
    if any(re.search(p, seq, re.IGNORECASE) for p in patterns):
        # trying to locate the repetitive domain (previously hvd) within the HYP-bearing reads
        hvd_all=find_hvd(before, after, seq, qual, match_thresh=0.25)
        # Then finally substituting the motifs 
        if hvd_all is not None:
            hvd=str(hvd_all[0])[18:len(str(hvd_all[0]))-18+1]
            protein=fuzzy_replace_short_strings(hvd, motif_filename, which="protein", errors=1)
            symbolic=fuzzy_replace_short_strings(hvd, motif_filename, which="name", errors=1)
            row={"id": name, "DNA": hvd_all[0], "protein": protein, "symbolic": symbolic, "phred":hvd_all[1], "meanphred":mean_read_qual(hvd_all[1])}
            parsedlist.append(row)

parsed=pd.DataFrame(parsedlist)

# writing a file of parsed HYP1 reads that includes DNA, amino acid, motif names, base qualities, etc.
parsed.to_csv("parsed_" + fastq_filename.split('/')[-1] + ".txt", sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar='\\')

