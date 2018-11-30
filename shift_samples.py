import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
parser.add_argument('-s', '--shift_length', type=int, action='store', dest='shift_length', help='Max number of bases to shift samples by')
args = parser.parse_args()

def shift_sample(seq,shift):
	if shift == 0:
		return seq
	else:
		return np.concatenate((seq[shift:,...],seq[:shift,...]))

# load sample metadata
L, Nc, ns, nb, nw, b = np.genfromtxt(f'../test_data/output_{args.infile}/{args.outfile}/samples_inf.txt',dtype='int32')
# Count samples
n_samples = Nc*ns//nw
# load samples
samples = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy',mode='r',dtype=int,shape=(n_samples,L,3))
# Make new directory
if not os.path.isdir(f'../test_data/output_{args.infile}/{args.outfile}_shifted_{args.shift_length}'):
	os.mkdir(f'../test_data/output_{args.infile}/{args.outfile}_shifted_{args.shift_length}')
# Create an array to write shifted samples into
samples_shifted = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}_shifted_{args.shift_length}/samples.npy',mode='w+',dtype=int,shape=(n_samples,L,3))
# Copy metadata
np.savetxt(f'../test_data/output_{args.infile}/{args.outfile}_shifted/samples_inf.txt',(L,Nc,ns,nb,nw,b))
# Generate random shifts
shifts = np.random.randint(-args.shift_length,args.shift_length+1,size=n_samples)
print(shifts.size)
# shift each sample by a random amount
for i in range(n_samples):
	samples_shifted[i,...] = shift_sample(samples[i,...],shifts[i])


