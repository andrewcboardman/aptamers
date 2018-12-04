import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-L', type=int, action='store',help='Number of spins')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
parser.add_argument('-s', type=int, action='store',default=20,help='Size of centre matrix')
args = parser.parse_args()

# Initiate J with random Gaussian distributed numbers, symmetrise, set diagonal to zero, and scale by 1/rt(N)
J = np.zeros((args.L,3,args.L,3))
J_cent = np.random.normal(size=(args.s,3,args.s,3))
gap = (args.L - args.s)//2
J[gap:(args.L-gap),:,gap:(args.L-gap),:] = J_cent

J = 0.25*(J + np.transpose(J,(2,1,0,3)) + np.transpose(J,(0,3,2,1)) + np.transpose(J,(2,3,0,1)))
J = J/np.sqrt(3*args.L)
J[range(args.L),:,range(args.L),:] = 0
h_cent = np.random.normal(size=(args.s,3))/np.sqrt((3*args.L))
h = np.concatenate((np.zeros((gap,3)),h_cent,np.zeros((gap,3))))

# Make output folder 
if not os.path.isdir('../test_data/output_{}'.format(args.outfile)):
	os.mkdir('../test_data/output_{}'.format(args.outfile))
# save 
np.save('../test_data/output_{}/fields_init.npy'.format(args.outfile),h)
np.save('../test_data/output_{}/couplings_init.npy'.format(args.outfile),J)

