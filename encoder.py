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


# Parse sequences and extract <s_i> and <s_i*s_j>

def averages(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# convert into a data vector
	mag1s = np.where(chars[:,np.newaxis] == ['A','C','G','T'],1,0).flatten()
	mag2s = np.outer(mag1s,mag1s)
	return (mag1s,mag2s)

mag1s_av = np.zeros(4*L)
mag2s_av = np.zeros((4*L,4*L))

seqs = SeqIO.parse(args.infile,'fasta')
mags = (averages(seq) for seq in seqs)
for mag in mags:
	mag1s_av += mag[0]/N
	mag2s_av += mag[1]/N

# correlations = matrix of dot products of columns - outer product of magnetisations
corr = mag2s_av - np.outer(mag1s_av,mag1s_av)*N

np.savetxt('magnetisations.out',mag1s_av)
np.savetxt('correlations.out',corr)




