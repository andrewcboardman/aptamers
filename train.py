
import numpy as np
import argparse
from Bio import SeqIO
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of flag for magnetisation and correlation files')
parser.add_argument('-m','--model-type', type=str, action='store', dest='model_type',default='mean_field', help='Type of model to be trained')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()


def parse_infile(infile):
	# Read fasta files
	if args.infile.endswith('fasta'):
		print('Reading {}...'.format(args.infile))
		return SeqIO.parse(args.infile,'fasta')
	else:
		print('Warning: unrecognised input file format')

def encode(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# One-hot encoding to find magnetisations
	return np.where(chars.reshape(L,1) == ['C','G','T'],1,0)

# Set up folder for output
if not os.path.isdir('../test_data/output_{}'.format(args.outfile)):
	os.mkdir('../test_data/output_{}'.format(args.outfile))

# Begin parsing input
input_parser = parse_infile(args.infile)

# Choose which training to use
if args.model_type == 'ind_site':
	# Independent site model- no couplings between sites
	# Parse input
	f1s = np.zeros((L,3))
	for seq in input_parser:
		f1s += encode(seq)/N
	# Add pseudocount
	f1s += 1/N
	# Fields are log(frequency) + a constant
	h_ind = np.log(f1s)
	np.save('../test_data/output_{}/h_ind.npy',h_ind)

elif args.model_type == 'mean_field':
	# Mean-field model
	# Factorise in sites - self-consistent
	# Parse input	
	f1s = np.zeros((L,3))
	f2s = np.zeros((L*3,L*3))
	print('Counting frequencies...')
	fs = ((seq,np.outer(seq,seq)) for seq in input_parser)
	for f in fs:
		f1s += f[0]/N
		f2s += f[1]/N

	# Add pseudocounts
	f1s += 1/N
	f2s += 1/N

	# correlations = average two-site magnetisations - outer product of one-site magnetisations
	corr = f2s - np.outer(f1s,f1s)
	print('Calculated frequencies...')

	# Use mean-field approximation to find fields and couplings
	print('Calculating fields and couplings')
	import matplotlib.pyplot as plt

	# Couplings inferred from inverse correlations
	J_mf = -np.linalg.inv(corr).reshape(L,3,L,3)

	# only off-diagonal couplings are important
	for i in range(L):
		J_mf[i,:,i,:] = 0

	# Infer fields
	h_mf = np.arctanh(f1s) - np.tensordot(J_mf,f1s,axes=2)

	# Save fields and couplings
	np.save('../test_data/output_{}/fields.npy'.format(args.outfile),h_mf)
	np.save('../test_data/output_{}/couplings.npy'.format(args.outfile),J_mf)

elif args.model_type == 'ising_ace':
	print("This hasn't been implemented yet lol")


