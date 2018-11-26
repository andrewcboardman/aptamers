import numpy as np
from Bio import SeqIO
import argparse
import itertools

# Take name of read file as input at command line

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of raw reads file (FASTQ format)')
parser.add_argument('-k', '--kmers', type=int, action='store',dest='k',default='1',help='k-mer size to use as feature')
args = parser.parse_args()
k = args.k


# Sequences must all be the same length and not contain any ambiguous base calls

seqs = SeqIO.parse(args.infile,'fasta')
N = 0
L = 0
for seq in seqs:
	N += 1
	L = max(L,len(seq))
print('{} reads'.format(N))
print('Length {}'.format(L))

# Generate all k-mers of length k present in the sequences

seqs = SeqIO.parse(args.infile,'fasta')
def KMerSplit(seq,k):
	return (seq[i:i+k] for i in range(L-k))
kmers = np.empty(N*(L-k),dtype=str)
i = 0
for seq in seqs:
	for kmer in KMerSplit(seq,k):
		kmers[i] = str(kmer)

# Find unique kmers
kmers_unique,kmers_index,kmer_counts = np.unique(kmers,return_index=True,return_counts=True)

print(np.max(kmers_counts))
helpI

bases = ('A','C','G','T')
kmers = (''.join(p) for p in itertools.product(bases, repeat=args.k))


kmerMag = np.zeros(4**args.k)

for i,kmer in enumerate(kmers):
	seqs = SeqIO.parse(args.infile,'fasta')
	for seq in seqs:
		kmerMag[i] += (kmer in seq)/N

print(kmerMag)