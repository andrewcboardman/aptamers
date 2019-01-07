import matplotlib.pyplot as plt
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='File containing energies of training sequences')
parser.add_argument('-c', '--control', type=str, action='store', dest='control', help='Flag containing energies of random sequences')
args = parser.parse_args()

energies = np.genfromtxt(args.path+args.infile)
plt.hist(energies,density=True,color='blue',alpha=0.5)
energies_rand = np.genfromtxt(args.path+args.control)
plt.hist(energies_rand,density=True,color='red',alpha=0.5)
plt.xlabel('Ising energy of sequence')
plt.ylabel('Frequency')
plt.show()
