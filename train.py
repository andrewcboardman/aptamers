
import numpy as np
import argparse
import os
import encode
import itertools as it
import load_samples

def IndSites(samples,n_spins,N):
	# Independent site model- no couplings between sites
	# Parse input
	f1s = np.zeros((n_spins,3))
	for sample in samples:
		f1s += sample/N
	# Add pseudocount
	f1s += 1/N
	# Fields are log(frequency) + a constant
	return np.log(f1s)

def MeanField(samples,n_spins,n_samples,n_states,correct):
	f1s = np.zeros((n_spins,n_states))
	f2s = np.zeros((n_spins*n_states,n_spins*n_states))
	print('Counting frequencies...')
	fs = ((seq,np.outer(seq,seq)) for seq in samples)
	for f in fs:
		f1s += f[0]/n_samples
		f2s += f[1]/n_samples
	# pseudocount to prevent overflow
	f1s = np.clip(f1s,None,0.999)
	# correlations = average two-site magnetisations - outer product of one-site magnetisations
	corr = f2s - np.outer(f1s,f1s)
	print('Calculated frequencies...')

	# Use mean-field approximation to find fields and couplings
	print('Calculating fields and couplings')
	import matplotlib.pyplot as plt

	# Couplings inferred from inverse correlations
	J_mf = -np.linalg.inv(corr).reshape(n_spins,n_states,n_spins,n_states)

	# only off-diagonal couplings are important
	for i in range(n_spins):
		J_mf[i,:,i,:] = 0

	# Infer fields using first-order approximation
	h_mf = np.log(f1s / np.clip((1 - np.sum(f1s,axis=1)),0.001,None).reshape(n_spins,1))
	if correct:
		h_mf -= np.tensordot(J_mf,f1s,axes=2)
	return (h_mf,J_mf)

def MeanFieldSlice(samples,n_spins,N,correct,slice_length):
	f1s = np.zeros((slice_length,3))
	f2s = np.zeros((slice_length*3,slice_length*3))
	print('Counting frequencies...')
	
	# use Cartesian product to iterate over all slices of the samples
	sliced_samples = (sample[i:(i+slice_length)] for (i,sample) in it.product(range(n_spins-slice_length),samples))
	fs = ((seq,np.outer(seq,seq)) for seq in sliced_samples)

	# Count one- and two-site frequencies
	for f in fs:
		f1s += f[0]/(N*(n_spins-slice_length))
		f2s += f[1]/(N*(n_spins-slice_length))

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
	parser.add_argument('-en', '--encoding', type=str,action='store',dest='encoding', help='Method of encoding sequence data')
	parser.add_argument('-s', '--slice_length', type=int, action='store', dest='slice_length', help='Length of slice to take of sample')
	parser.add_argument('-k', type=int, action='store',help='length of kmer used')
	args = parser.parse_args()

	# load sample metadata
	n_samples,n_spins,b = load_samples.load_metadata(f'../test_data/output_{args.infile}/{args.outfile}/samples_inf.txt')

	# load samples
	if args.encoding == 'bases':
		n_states = 3
	elif args.encoding == 'bases_align':
		n_states = 4
	else:
		n_states = 1
		n_spins = 4**args.k

	samples = load_samples.load_samples(f'../test_data/output_{args.infile}/{args.outfile}/samples_{args.encoding}.txt',n_samples,n_spins,n_states)

	# execute training routine
	if args.model == 'mean_field':
		(h_mf,J_mf) = MeanField(samples,n_spins,n_samples,n_states,args.correct) # scale by inverse T
		h_mf = h_mf/b
		J_mf = J_mf/b
		# Save fields and couplings
		np.savetxt(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf.txt',h_mf)
		np.savetxt(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf.txt',J_mf.reshape(n_spins*n_states,n_spins*n_states))
	elif args.model == 'ind_sites':
		h_ind = IndSites(samples,n_samples)
		np.save(f'../test_data/output_{args.outfile}/field_ind.npy',h_ind)
	elif args.model == 'mean_field_slice':
		(h_mf,J_mf) = MeanFieldSlice(samples,n_spins,Nc*ns//nw,args.correct,args.slice_length) # scale by inverse T
		h_mf = h_mf/b 
		J_mf = J_mf/b
		# Save fields and couplings
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf_slice_{args.slice_length}.npy',h_mf)
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf_slice_{args.slice_length}.npy',J_mf)


if __name__ == '__main__':
	main()