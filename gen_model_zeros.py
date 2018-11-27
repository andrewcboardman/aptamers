import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-L', type=int, action='store',help='Number of spins')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()

# Initiate J with random Gaussian distributed numbers, symmetrise, set diagonal to zero, and scale by 1/rt(N)
J = np.zeros((args.L,3,args.L,3))
J = 0.25*(J + np.transpose(J,(2,1,0,3)) + np.transpose(J,(0,3,2,1)) + np.transpose(J,(2,3,0,1)))
J = J/np.sqrt(3*args.L)
J[range(args.L),:,range(args.L),:] = 0
h = np.zeros((args.L,3))

# Make output folder 
if not os.path.isdir('../test_data/output_{}'.format(args.outfile)):
	os.mkdir('../test_data/output_{}'.format(args.outfile))
# save 
np.save('../test_data/output_{}/fields_init.npy'.format(args.outfile),h)
np.save('../test_data/output_{}/couplings_init.npy'.format(args.outfile),J)

