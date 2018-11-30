import numpy as np
import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def seq2array(seq,L):
	array = np.zeros((L,3))
	seqarr = np.array(list(seq))[:,np.newaxis]
	array[seqarr==np.array(('C','G','T'))] = 1
	return array

def encode(seqs,arr,L):
	strings = (str(seq.seq) for seq in seqs)
	for i,string in enumerate(strings):
		arr[i,...] = seq2array(string,L)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
	args = parser.parse_args()

	# Read FASTA-formatted samples
	seqs = SeqIO.parse(f'../test_data/output_{args.infile}/{args.outfile}/samples.fasta','fasta')

	# load sample metadata
	L, Nc, ns, nb, nw, b = np.genfromtxt(f'../test_data/output_{args.infile}/{args.outfile}/samples_inf.txt',dtype='int32')

	# Count samples
	n_samples = Nc*ns//nw

	# Make output array for samples
	samples = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy',mode='w+',dtype=int,shape=(n_samples,L,3))


	# Populate array with encoded strings
	encode(seqs,samples,L)

	# Close connection to array
	del samples

if __name__ == '__main__':
	main()