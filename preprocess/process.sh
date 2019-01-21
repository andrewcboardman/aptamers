#!/usr/bin/bash

# Select reads containing forward adapters + length 40 sequence
cutadapt minus_IgE_same_PE_R2.fastq -g AGCAGCACAGAGGTCAGATG --discard-untrimmed -o minusR2_trim_forward_tmp.fastq
cutadapt minusR2_trim_forward_tmp.fastq -a CCTATGCGTGCTACCGTGAA --discard-untrimmed -o minusR2_trim_forward_tmp2.fastq
cutadapt minusR2_trim_forward_tmp2.fastq -m 40 -M 40 -o minusR2_trim_forward.fastq
# Select reads containing backward adapters + length 40 sequence
cutadapt minus_IgE_same_PE_R2.fastq -g TTCACGGTAGCACGCATAGG --discard-untrimmed -o minusR2_trim_backward_tmp.fastq
cutadapt minusR2_trim_backward_tmp.fastq -a CATCTGACCTCTGTGCTGCT --discard-untrimmed -o minusR2_trim_backward_tmp2.fastq
cutadapt minusR2_trim_backward_tmp2.fastq -m 40 -M 40 -o minusR2_trim_backward.fastq

paste - - - - < plusR2_trim_forward.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > plusR2_trim_forward.fasta
paste - - - - < plusR2_trim_backward.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > plusR2_trim_backward.fasta

fastx_reverse_complement -i plusR2_trim_backward.fasta -o plusR2_trim_backward_rev.fasta

cat minusR2_trim_forward.fasta minusR2_trim_backward_rev.fasta > minusR2_combine.fasta

fastx_collapser -i minusR2_combine.fasta -o minusR2.fasta

cat minusR2.fasta minusR1.fasta > minus_combine.fasta

fastx_collapser -i minus_combine.fasta -o minus.fasta