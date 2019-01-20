import matplotlib.pyplot as plt
import argparse
import numpy as np

def plot(E1,E2):
	plt.hist(E1,density=True,color='blue',alpha=0.5)
	plt.hist(energies_rand,density=True,color='red',alpha=0.5)
	plt.xlabel('Ising energy of sequence')
	plt.ylabel('Frequency')
	plt.show()



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')
	parser.add_argument('-i1', '--infile1', type=str, action='store', dest='infile1',	help='File containing energies of first set of sequences')
	parser.add_argument('-i2', '--infile2', type=str, action='store', dest='infile2',	help='File containing energies of second set of sequences')
	args = parser.parse_args()
	energies = np.genfromtxt(args.path+args.infile1)
	energies_rand = np.genfromtxt(args.path+args.infile2)
	plot(energies,energies_rand)

if __name__ = '__main__':
	main()