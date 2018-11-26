import numpy as np 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of samples file')
parser.add_argument('-c', '--cutoff', type=str, action='store', dest='cutoff',default='250000',	help='number of samples to drop from beginning of file')
args = parser.parse_args()
samples = np.genfromtxt('{}.out'.format(args.infile))
np.save('{}_trimmed.npy'.format(args.infile),samples[cutoff:,:])

