import argparse
import numpy as np
from scipy.signal import correlate
import matplotlib.pyplot as plt

def PlotAutoCorr(samples,n,L):
	# normalise samples
	normed_samples = (samples - np.mean(samples,axis=0))/np.std(samples,axis=0)
	# calculate autocorrelation
	corr = correlate(normed_samples,normed_samples,mode='same',method='fft')/(n*L*3)
	return plt.plot(range(-n//2,n//2),corr[:,L//2,1])


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',help='Flag for original model')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
	parser.add_argument('-a','--all',action='store_true',help='Calculate all autocorrelations (rather than just 1)')
	args = parser.parse_args()

	
	# load options
	#L, Nc, ns, nb, nw, b = np.genfromtxt(f'../test_data/output_{args.infile}/{args.outfile}/samples_inf.txt',dtype='int32')
	#n_samples = ns//nw
	# parse npy file using mmap loading so we don't overload RAM
	#samples = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy',mode='r',dtype=int,shape=(n_samples,Nc,L,3))
	if args.all:
		for i in range(Nc):
			PlotAutoCorr(samples[:,i,...],n_samples,L)
	else:
		PlotAutoCorr(samples[:,np.random.randint(Nc),...],n_samples,L)

	plt.xlabel(f'Offset / {nw} MC steps')
	plt.ylabel('Normalised autocorrelation')
	plt.show()

if __name__ == '__main__':
	main()