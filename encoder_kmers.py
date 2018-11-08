import numpy as np
from Bio import SeqIO
import argparse
import itertools

# Take name of read file as input at command line

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of raw reads file (FASTQ format)')
parser.add_argument('-k', '--kmers', type=int, action='store',dest='k',default='1',help='k-mer size to use as feature')
args = parser.parse_args()

# Sequences must all be the same length and not contain any ambiguous base calls

seqs = SeqIO.parse(args.infile,'fasta')
N = 0
L = 0
for seq in seqs:
	N += 1
	L = max(L,len(seq))
print('{} reads'.format(N))
print('Length {}'.format(L))

# Generate all k-mers of length k 

bases = ('A','C','G','T')
kmers = (''.join(p) for p in itertools.product(bases, repeat=args.k))


kmerMag = np.zeros(4**args.k)

for i,kmer in enumerate(kmers):
	seqs = SeqIO.parse(args.infile,'fasta')
	for seq in seqs:
		kmerMag[i] += (kmer in seq)/N

print(kmerMag)