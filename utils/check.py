import argparse
import numpy as np
import load_samples


def main():
	"""Calculate the probability of a set of sequences given a set of coefficients"""
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Samples filename')
	parser.add_argument('-m','--model', type=str, action='store', dest='model', help='model coefficients file tag')
	parser.add_argument('-en', '--encoding', type=str,action='store',dest='encoding', help='Method of encoding sequence data')

	parser.add_argument('-N', '--n_samples', type=int, action='store', help='Number of samples')
	parser.add_argument('-q', '--n_states', type=int, action='store', help='Number of states per spin')
	parser.add_argument('-L', '--n_spins', type=int, action='store', help='Number of spins')
	parser.add_argument('-k', type=int, action='store',help='length of kmer used')
	args = parser.parse_args()


	samples = load_samples.load_samples(args.path+args.infile,args.n_samples,args.n_spins,args.n_states)

	# load model
	h_mf = np.genfromtxt(args.path + 'fields_mf.txt')
	J_mf = np.genfromtxt(args.path + 'couplings_mf.txt').reshape(args.n_spins,args.n_states,args.n_spins,args.n_states)

	# calculate energies
	energies = np.zeros(args.n_samples)
	energies -= np.tensordot(samples,h_mf,axes=2)
	energies -= np.einsum('ijk,jklm,ilm->i',samples,J_mf,samples)/2

	np.savetxt(args.path+args.infile+args.model,energies)


if __name__ == '__main__':
	main()