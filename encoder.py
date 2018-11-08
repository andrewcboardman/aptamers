import numpy as np
from Bio import SeqIO
import argparse

# Take name of read file as input at command line

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of raw reads file (FASTQ format)')
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


# Create encoded array

def encode(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# convert into a data vector
	boole = np.where(chars[:,np.newaxis] == ['A','C','G'],1,-1).flatten()
	return boole

seqs = SeqIO.parse(args.infile,'fasta')
coded_array = np.zeros((N,3*L))
for i,rec in enumerate(seqs):
	coded_array[i,:] = encode(rec)

# Calculate the magnetisations and connected correlations of the samples
# magnetisations = average column value
mag = np.sum(coded_array,axis=1)/N
# correlations = matrix of dot products of columns - outer product of magnetisations
corr = np.tensordot(coded_array,coded_array,axes=(1,1))/N - np.outer(mag,mag)

np.savetxt('magnetisations.out',mag)
np.savetxt('correlations.out',corr)




