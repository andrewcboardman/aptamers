import argparse
import numpy as np
import load_samples

def main():
	"""Calculate the probability of a set of sequences given a set of coefficients"""
	parser = argparse.ArgumentParser()
	parser.add_argument('-m','--model', type=str, action='store', dest='model',default='mean_field', help='Type of model to be trained')
	parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Samples')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
	parser.add_argument('-en', '--encoding', type=str,action='store',dest='encoding', help='Method of encoding sequence data')
	parser.add_argument('-k', type=int, action='store',help='length of kmer used')
	args = parser.parse_args()

	# load sample metadata
	n_samples,n_spins,b = load_samples.load_metadata(args.path + 'samples_inf.txt')
	# load samples
	if args.encoding == 'bases':
		n_states = 3
	elif args.encoding == 'bases_align':
		n_states = 4
	else:
		n_states = 1
		n_spins = 4**args.k

	samples = load_samples.load_samples(args.path+args.infile,n_samples,n_spins,n_states)

	# load model
	h_mf = np.genfromtxt(args.path + 'fields_mf.txt')
	J_mf = np.genfromtxt(args.path + 'couplings_mf.txt').reshape(n_spins,n_states,n_spins,n_states)

	# calculate energies
	energies = np.zeros(n_samples)
	energies -= np.tensordot(samples,h_mf,axes=2)
	energies -= np.einsum('ijk,jklm,ilm->i',samples,J_mf,samples)/2

	np.savetxt(args.path+args.outfile,energies)


if __name__ == '__main__':
	main()