
import numpy as np
import argparse
from Bio import SeqIO
import os
import encode
import itertools as it

def IndSites(samples,L,N):
	# Independent site model- no couplings between sites
	# Parse input
	f1s = np.zeros((L,3))
	for sample in samples:
		f1s += sample/N
	# Add pseudocount
	f1s += 1/N
	# Fields are log(frequency) + a constant
	return np.log(f1s)

def MeanField(samples,L,N,correct):
	f1s = np.zeros((L,3))
	f2s = np.zeros((L*3,L*3))
	print('Counting frequencies...')
	fs = ((seq,np.outer(seq,seq)) for seq in samples)
	for f in fs:
		f1s += f[0]/N
		f2s += f[1]/N
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

	# Infer fields using first-order approximation
	h_mf = np.log(f1s / (1 - np.sum(f1s,axis=1)).reshape(L,1))
	if correct:
		h_mf -= np.tensordot(J_mf,f1s,axes=2)
	return (h_mf,J_mf)

def MeanFieldSlice(samples,L,N,correct,slice_length):
	f1s = np.zeros((slice_length,3))
	f2s = np.zeros((slice_length*3,slice_length*3))
	print('Counting frequencies...')
	
	# use Cartesian product to iterate over all slices of the samples
	sliced_samples = (sample[i:(i+slice_length)] for (i,sample) in it.product(range(L-slice_length),samples))
	fs = ((seq,np.outer(seq,seq)) for seq in sliced_samples)

	# Count one- and two-site frequencies
	for f in fs:
		f1s += f[0]/(N*(L-slice_length))
		f2s += f[1]/(N*(L-slice_length))

	# correlations = average two-site magnetisations - outer product of one-site magnetisations
	corr = f2s - np.outer(f1s,f1s)
	print('Calculated frequencies...')

	# Use mean-field approximation to find fields and couplings
	print('Calculating fields and couplings')
	import matplotlib.pyplot as plt

	# Couplings inferred from inverse correlations
	J_mf = -np.linalg.inv(corr).reshape(slice_length,3,slice_length,3)

	# only off-diagonal couplings are important
	for i in range(slice_length):
		J_mf[i,:,i,:] = 0

	# Infer fields using first-order approximation
	h_mf = np.log(f1s / (1 - np.sum(f1s,axis=1)).reshape(slice_length,1))
	if correct:
		h_mf -= np.tensordot(J_mf,f1s,axes=2)
	return (h_mf,J_mf)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m','--model', type=str, action='store', dest='model',default='mean_field', help='Type of model to be trained')
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
	parser.add_argument('-c', '--correct', action='store_true', dest='correct', help='Apply correction to fields')
	parser.add_argument('-s', '--slice_length', type=int, action='store', dest='slice_length', help='Length of slice to take of sample')
	parser.add_argument('-r', '--remove', action='store_true', dest='remove', help='remove numpy binary')
	args = parser.parse_args()

	# load sample metadata
	metadata = np.genfromtxt(f'../test_data/output_{args.infile}/{args.outfile}/samples_inf.txt',dtype='float32')
	L, Nc, ns, nb, nw = metadata[:-1].astype('int32')
	b = metadata[-1]

	# Read FASTA-formatted samples
	#seqs = SeqIO.parse(f'../test_data/output_{args.infile}/{args.outfile}/samples.fasta','fasta')
	
	# Create numpy array (stored on hard drive to prevent RAM overloading) 
	samples = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy',mode='r',dtype=int,shape=(Nc*ns//nw,L,3))

	# Populate numpy array with encoded samples
	#encode.encode(seqs,samples,L)

	# execute training routine
	if args.model == 'mean_field':
		(h_mf,J_mf) = MeanField(samples,L,Nc*ns//nw,args.correct) # scale by inverse T
		h_mf = h_mf/b 
		J_mf = J_mf/b
		# Save fields and couplings
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf.npy',h_mf)
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf.npy',J_mf)
	elif args.model == 'ind_sites':
		h_ind = IndSites(samples,L,Nc*ns//nw)
		np.save(f'../test_data/output_{args.outfile}/field_ind.npy',h_ind)
	elif args.model == 'mean_field_slice':
		(h_mf,J_mf) = MeanFieldSlice(samples,L,Nc*ns//nw,args.correct,args.slice_length) # scale by inverse T
		h_mf = h_mf/b 
		J_mf = J_mf/b
		# Save fields and couplings
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf_slice_{args.slice_length}.npy',h_mf)
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf_slice_{args.slice_length}.npy',J_mf)

	del samples
	if args.remove:
		os.remove(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy')

if __name__ == '__main__':
	main()