import numpy as np
from Bio import SeqIO
import argparse

# Take name of read file as input at command line

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, action='store', dest='infile',default='test.fastq'
	help='Name of raw reads file (FASTQ format)')
args = parser.parse_args()


# Find max read length and number of reads

n=0
l=0
for rec in SeqIO.parse(infile,'fastq'):
	n +=1
	if len(str(rec)) > l:
		l = len(str(rec))

# Create array and fill with encoding

def encode(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# convert into a data vector
	boole = 1*(chars[:,np.newaxis] == ['A','C','G']).flatten()
	# pad to get to vector with length 120
	pad_boole = np.pad(boole,(0,120-len(boole)))
	# Join vector to string
	return pad_boole

coded_array = np.zeros((n,120))
for i,rec in enumerate(SeqIO.parse(infile,'fastq')):
	coded_array[i,:] = encode(rec)

# Calculate magnetisations and connected correlations
m = np.sum(coded_array,axis=1)/n
s = np.tensordot(coded_array,coded_array,axes=1)/n - np.outer(m,m)

np.savetext('magnet.out',m)
np.savetext('correlations.out',s)




