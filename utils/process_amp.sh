#!/usr/bin/bash
FILE=$1

fastx_quality_stats -i ${FILE}.fastq -o ${FILE}.init.stats #nucleotide distribution and quality scores


cutadapt -v -a AGCAGCACAGAGGTCAGATG...CCTATGCGTGCTACCGTGAA ${FILE}.forward.fastq | \
fastq_quality_filter -v -q 20 -p 100 | \
fastq_to_fasta | \
lenfilter.py -L 40 | \
fastx_collapser -v -o ${FILE}_process.forward.fasta

cutadapt -v -a TTCACGGTAGCACGCATAGG...CATCTGACCTCTGTGCTGCT ${FILE}.reverse.fastq | \
fastq_quality_filter -v -q 20 -p 100 | \
fastq_to_fasta | \
lenfilter.py -L 40 | \
fastx_collapser -v -o ${FILE}_process.reverse.fasta


fastx_quality_stats -i ${FILE}_process.forward.fasta -o ${FILE}.post.forward.stats
fastx_quality_stats -i ${FILE}_process.reverse.fasta -o ${FILE}.post.reverse.stats