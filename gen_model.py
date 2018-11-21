import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()


J_mf = np.zeros((120,120))
for i in range(0,120,3):
	for j in range(0,120,3):
		if i !=
		J_mf[i:i+3,j:j+3] =np.identity(3)

h_mf = np.ones(120)
import matplotlib.pyplot as plt
plt.imshow(J_mf)
plt.show()
if not os.path.isdir('../test_data/output_{}'.format(args.outfile)):
	os.mkdir('../test_data/output_{}'.format(args.outfile))

np.savetxt('../test_data/output_{}/fields_{}.out'.format(args.outfile,args.outfile),h_mf)
np.savetxt('../test_data/output_{}/couplings_{}.out'.format(args.outfile,args.outfile),J_mf)

