
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

def MeanField(f1s,f2s,n_spins,n_states,correct):
	
	# pseudocount to prevent overflow
	f1s = np.clip(f1s,None,0.999)
	# correlations = average two-site magnetisations - outer product of one-site magnetisations
	corr = f2s - np.outer(f1s,f1s)

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
	
	parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='File containing encoded samples')
	
	parser.add_argument('-c', '--correct', action='store_true', dest='correct', help='Apply correction to fields')
	
	parser.add_argument('-q', '--n_states', type=str, action='store', help='Number of states per spin')
	parser.add_argument('-L', '--n_spins', type=str, action='store', help='Number of spins')
	parser.add_argument('-b', '--temp', type=str, action='store', help='Inverse temperature')
	args = parser.parse_args()

	f1s = np.genfromtxt(args.path + 'f1s_' + args.infile)
	f2s = np.genfromtxt(args.path + 'f2s_' + args.infile)

	# execute training routine
	if args.model == 'mean_field':
		(h_mf,J_mf) = MeanField(f1s,f2s,n_spins,n_states,args.correct) # scale by inverse T
		h_mf = h_mf/b
		J_mf = J_mf/b
		# Save fields and couplings
		np.savetxt(args.path + 'h_mf_' + args.infile,h_mf)
		np.savetxt(args.path + 'J_mf_' + args.infile,J_mf.reshape(n_spins*n_states,n_spins*n_states))
	elif args.model == 'ind_sites':
		h_ind = IndSites(samples,n_samples)
		np.savetxt(args.path + 'fields_ind.txt',h_ind)
	elif args.model == 'mean_field_slice':
		(h_mf,J_mf) = MeanFieldSlice(samples,n_spins,Nc*ns//nw,args.correct,args.slice_length) # scale by inverse T
		h_mf = h_mf/b 
		J_mf = J_mf/b
		# Save fields and couplings
		np.savetxt(args.path + 'fields_mf_slice.txt',h_mf)
		np.savetxt(args.path + 'couplings_mf_slice.txt',J_mf.reshape(args.slice_length*n_states,args.slice_length*n_states))


if __name__ == '__main__':
	main()