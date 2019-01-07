import matplotlib.pyplot as plt
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
args = parser.parse_args()
path = f'../test_data/output_{args.infile}/{args.outfile}/'

energies = np.genfromtxt(path+'couplings_mf.txt').flatten()
plt.hist(energies,density=True,bins=50)
plt.xlabel('Coupling strength')
plt.ylabel('Frequency')
plt.show()
