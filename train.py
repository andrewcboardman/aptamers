
import numpy as np
import argparse
from Bio import SeqIO
import os

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

def MeanField(samples,L,N):
	f1s = np.zeros((L,3))
	f2s = np.zeros((L*3,L*3))
	print('Counting frequencies...')
	fs = ((seq,np.outer(seq,seq)) for seq in samples)
	for f in fs:
		f1s += f[0]/N
		f2s += f[1]/N
	print(np.mean(np.sum(f1s,axis=1)-0.75))
	print(f2s)
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
	h_mf = np.log(f1s / (1 - np.sum(f1s,axis=1)).reshape(L,1)) - np.tensordot(J_mf,f1s,axes=2)
	print(np.tensordot(J_mf,f1s,axes=2))

	return (h_mf,J_mf)



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m','--model', type=str, action='store', dest='model',default='mean_field', help='Type of model to be trained')
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
	args = parser.parse_args()

	# load sample metadata
	L, Nc, ns, nb, nw, b = np.genfromtxt(f'../test_data/output_{args.infile}/{args.outfile}/samples_inf.txt',dtype='int32')

	# parse npy file using mmap loading so we don't overload RAM
	# merge first two columns 
	samples = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy',mode='r',dtype=int,shape=(Nc*ns//nw,L,3))
	#samples = np.load(f'../test_data/output_{args.outfile}/samples.npy',mmap_mode='r').reshape(Nc*ns//nw,L,3)
	# execute training routine
	if args.model == 'mean_field':
		(h_mf,J_mf) = MeanField(samples,L,Nc*ns//nw)
		# Save fields and couplings
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf.npy',h_mf)
		np.save(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf.npy',J_mf)
	elif args.model == 'ind_sites':
		h_ind = IndSites(samples,L,Nc*ns//nw)
		np.save(f'../test_data/output_{args.outfile}/field_ind.npy',h_ind)


if __name__ == '__main__':
	main()






