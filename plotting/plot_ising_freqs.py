import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Folder containing original couplings')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
parser.add_argument('-s', '--slice_length',type=int, action='store', dest='slice_length',default=0, help='length of slice from sequence used in training')
args = parser.parse_args()

def freq(h,J,L):
	f_start = np.exp( - h)
	f1 = np.exp( - h - 0.5 * np.tensordot(J, f_start, axes=2))
	f2 = np.exp( - h.reshape(L,3,1,1) - h.reshape(1,1,L,3) - 0.5 * J \
		- 0.5 * np.tensordot(J, f_start, axes=2).reshape(L,3,1,1) \
		- 0.5 * np.tensordot(J, f_start, axes=2).reshape(1,1,L,3) )
	return (f1,f2)


J_init = np.load(f'../test_data/output_{args.infile}/couplings_init.npy')
h_init = np.load(f'../test_data/output_{args.infile}/fields_init.npy')
L = h_init.shape[0]
f_init_1, f_init_2 = freq(h_init,J_init,L)


if args.slice_length == 0:
	J_mf = np.load(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf.npy')
	h_mf = np.load(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf.npy')
	f_mf_1, f_mf_2 = freq(h_mf,J_mf,L)

	fig, (ax1,ax2) = plt.subplots(1,2)

	ax1.scatter(f_init_1.flatten(),f_mf_1.flatten(),color='red',alpha=0.5)
	ax1.legend()
	ax1.set_xlabel(r'f$_{1s init}$')
	ax1.set_ylabel(r'f$_{1s MF}$')


	ax2.scatter(f_init_2.flatten(),f_mf_2.flatten(),color='red',alpha=0.5)
	ax2.legend()
	ax2.set_xlabel(r'f$_{2s init}$')
	ax2.set_ylabel(r'f$_{2s MF}$')

	plt.suptitle('Actual vs. inferred values for mean-field inference of a Gaussian matrix padded with zeros')
	plt.show()
else:
	J_mf_slice = np.load(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf_slice_{args.slice_length}.npy')
	h_mf_slice = np.load(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf_slice_{args.slice_length}.npy')

	L = J_init.shape[0]
	gap = (L - args.slice_length) //2
	J_mf = np.zeros((L,3,L,3))
	J_mf[gap:(L-gap),:,gap:(L-gap),:] = J_mf_slice
	h_mf = np.zeros((L,3))
	h_mf[gap:(L-gap),:] = h_mf_slice

	f_mf_1, f_mf_2 = freq(h_mf,J_mf,L)

	fig, (ax1,ax2) = plt.subplots(1,2)

	ax1.scatter(f_init_1.flatten(),f_mf_1.flatten(),color='red',alpha=0.5)
	ax1.legend()
	ax1.set_xlabel(r'f$_{1s init}$')
	ax1.set_ylabel(r'f$_{1s MF}$')


	ax2.scatter(f_init_2.flatten(),f_mf_2.flatten(),color='red',alpha=0.5)
	ax2.legend()
	ax2.set_xlabel(r'f$_{2s init}$')
	ax2.set_ylabel(r'f$_{2s MF}$')


	plt.suptitle('Actual vs. inferred values for mean-field inference of a Gaussian matrix padded with zeros')
	plt.show()