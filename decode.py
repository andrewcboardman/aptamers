import numpy as np
import argparse
import os
from Bio.SeqIO import write
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def decode(sample,L):
	seq = np.empty(dtype='U1',shape=L)
	seq[sample[:,0]==1] = 'C'
	seq[sample[:,1]==1] = 'G'
	seq[sample[:,2]==1] = 'T'
	seq[seq==''] = 'A'
	return ''.join(seq)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
	args = parser.parse_args()

	# load sample metadata
	L, Nc, ns, nb, nw, b = np.genfromtxt(f'../test_data/output_{args.infile}/{args.outfile}/samples_inf.txt',dtype='int32')

	# Count samples
	n_samples = Nc*ns//nw

	# load samples
	samples = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy',mode='r',dtype=int,shape=(n_samples,L,3))

	# Decode to string
	strings = (decode(sample,L) for sample in samples)

	# Convert to SeqRecord for writing to file
	seqs = (SeqRecord(Seq(string,IUPAC.IUPACUnambiguousDNA()),id=f'sample_{i}') for (i,string) in enumerate(strings))

	# Write to FASTA
	write(seqs,f'../test_data/output_{args.infile}/{args.outfile}/samples.fasta','fasta')

if __name__ == '__main__':
	main()