import argparse
import numpy as np
from scipy.signal import correlate
import matplotlib.pyplot as plt

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of flag for magnetisation and correlation files')
	parser.add_argument('-m','--model', type=str, action='store', dest='model',default='mean_field', help='Type of model to be trained')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
	args = parser.parse_args()

	
	
	# load options
	L, Nc, ns, nb, nw, b = np.genfromtxt(f'../test_data/output_{args.outfile}/samples_10m_inf.txt',dtype='int32')
	n_samples = ns//nw
	# parse npy file using mmap loading so we don't overload RAM
	samples = np.memmap(f'../test_data/output_{args.outfile}/samples_10m.npy',mode='r',shape=(n_samples,Nc,L,3))
	


	for i in range(int(Nc)):
		mean_samples = np.mean(samples[:,i,...],axis = 0)
		std_samples = np.std(samples[:,i,...],axis = 0)
		normed_samples = (samples[:,i,...] - mean_samples)
		corr = correlate(normed_samples,normed_samples,mode='full',method='fft')/(ns//nw)
		plt.scatter(range(ns//nw-1),corr[ns//nw:,19,2])
	plt.show()

if __name__ == '__main__':
	main()