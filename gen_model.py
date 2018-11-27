import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-L', type=int, action='store',help='Number of spins')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()

# Initiate J as random Gaussian matrix
J = np.random.normal(size=(3*args.L,3*args.L)) # random values
J = 0.5*(J+J.T)/np.sqrt(3*args.L) # Scaled and symmetrised
J = J.reshape(args.L,3,args.L,3)
h = np.zeros((args.L,3))

# Make output folder 
if not os.path.isdir('../test_data/output_{}'.format(args.outfile)):
	os.mkdir('../test_data/output_{}'.format(args.outfile))
# save 
np.save('../test_data/output_{}/fields_init.npy'.format(args.outfile),h)
np.save('../test_data/output_{}/couplings_init.npy'.format(args.outfile),J)

