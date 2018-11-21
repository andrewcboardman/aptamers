import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of flag for magnetisation and correlation files')
parser.add_argument('-m','--model-type', type=str, action='store', dest='model_type',default='ising_mf', help='Type of model to be trained')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()

if args.model_type == 'ising_mf':
	
	# Load magnetisations and correlations from file
	h_mf = np.genfromtxt('../test_data/output_{}/fields_{}.out'.format(args.infile,args.infile))
	J_mf = np.genfromtxt('../test_data/output_{}/couplings_{}.out'.format(args.infile,args.infile))
	for i in range(120):
		J_mf[i:(i+3),i:(i+3)] = 0
	(fig,axs) = plt.subplots(1,2)
	axs[0].plot(h_mf)
	im = axs[1].imshow(J_mf)
	fig.colorbar(im)
	plt.show()