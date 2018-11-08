import numpy as np
from Bio import SeqIO
import argparse

# Take name of read file as input at command line

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of raw reads file (FASTQ format)')
args = parser.parse_args()

# Sequences must all be the same length and not contain any ambiguous base calls

seqs = SeqIO.parse(args.infile,'fastq')
N = 0
for sequence in seqs:
	N += 1
print('{} reads'.format(N))
L = len(seqs[1].seq)

# Create encoded array

def match(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# Change bases for integers mod 4
	match = (chars[:,np.newaxis] == ['A','C','G','T','-']).flatten()
	return match

def bias(seq):
	# Calculate the site-specific magnetisations and pairwise products for a sequence
	match = match(seq)
	mag = match/N
	pair = np.outer(match,match)/N
	return mag,pair


mag = np.zeros(5*L)
pair = np.zeros(5*L,5*L)

for seq in seqs:
	bias = bias(seq)
	mag += bias[0]
	pair += bias[1]

corr = pair - np.outer(match,match)

# Calculate the magnetisations and connected correlations of the samples
# magnetisations = average column value
mag = np.sum(coded_array[...,np.newaxis]==range(4),axis=1)/N

# correlations = matrix of dot products of columns - outer product of magnetisations
corr = np.tensordot(coded_array,coded_array,axes=(1,1))/N - np.outer(mag,mag)

np.savetxt('magnetisations.out',mag)
np.savetxt('correlations.out',corr)
