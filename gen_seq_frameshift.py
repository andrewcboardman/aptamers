import numpy as np
import argparse
import time

def IsingEnergy(seqs,h,J,N,L):
	"""Finds the Ising energy for each sample sequence"""
	field_energy = -np.tensordot(h,seqs,axes=((0,1),(1,2)))
	coupling_energy = np.diag(-np.tensordot(seqs,np.tensordot(seqs,J,axes=((1,2),(2,3))),axes=((1,2),(1,2))))
	return field_energy + coupling_energy

def GibbsSample(seqs,h,J,N,L,b):
	# Propose changes to sequences
	which_pos = np.random.randint(L,size=N)
	old_bases = seqs[range(N),which_pos,:]
	new_bases = 1*(np.random.randint(4,size=(N,1))==np.arange(1,4))
	base_changes = new_bases - old_bases

	# Calculate changes in energy using Ising model
	field_energy_changes = -np.sum(h[which_pos,:] * base_changes, axis=1) 
	coupling_energy_changes = -2*np.sum(np.sum(J[which_pos,...] * base_changes[...,np.newaxis,np.newaxis], axis=1) * seqs, axis=(1,2)) 
	deltas = field_energy_changes + coupling_energy_changes

	# Accept or reject changes based on energy change
	accepted = np.exp(-deltas*b)>np.random.rand(N)
	accepted_no = np.arange(N)[accepted]
	accepted_pos = which_pos[accepted]
	accepted_changes = base_changes[accepted]
	new_seqs = np.copy(seqs)
	new_seqs[accepted_no,accepted_pos,:] += accepted_changes
	return new_seqs

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
	parser.add_argument('-L', type=int, action='store',help='Number of spins (each have 4 states)',default=20)
	parser.add_argument('-Nc', type=int, action='store',help='Number of Markov chains',default=1)
	parser.add_argument('-ns', type=int, action='store',help='Number of sampling steps per chain',default=10000)
	parser.add_argument('-nb', type=int, action='store',help='Number of burn-in steps',default=1000)
	parser.add_argument('-nw', type=int, action='store',help='Number of steps between samples',default=100)
	parser.add_argument('-sp', action='store_true',help='Show plot of energies to check convergence',default=False)
	parser.add_argument('-b', type=float,action='store',help='Scale of energies (inverse temperature)',default=1)
	args = parser.parse_args()

	t0 = time.time()

	# Load fields and couplings from file
	h = np.load('../test_data/output_{}/fields.npy'.format(args.outfile))
	J = np.load('../test_data/output_{}/couplings.npy'.format(args.outfile))

	# Generate random sequences to initialise
	seqs = 1*(np.random.randint(4,size=(args.Nc,args.L,1))==np.arange(1,4))

	# Initialise output arrays
	n_samples = args.ns//args.nw
	energies = np.zeros((n_samples,args.Nc))

	# Samples is a rather large array so will be stored on disc as it is written (using memmap)
	samples = np.memmap(f'../test_data/output_{args.outfile}/samples_fs.mmap',dtype=int,mode='w+',shape=(n_samples,args.Nc,args.L,3))

	t1 = time.time()
	# Burn in 
	print('Started warmup...')
	for j in range(args.nb):
		seqs = GibbsSample(seqs,h,J,args.Nc,args.L,args.b)
	print('Warmup over.')
	print('Started sampling...')
	t2 = time.time()

	# Sample from Markov chain
	for k in range(args.ns):
		if k%args.nw == 0:
			print('{0}/{1} steps complete'.format(k,args.ns))
			samples[k//args.nw,...] = seqs
			energies[k//args.nw,:] = IsingEnergy(seqs,h,J,args.Nc,args.L)
		seqs = GibbsSample(seqs,h,J,args.Nc,args.L,args.b)
	t3 = time.time()

	# Save samples and energies
	np.save(f'../test_data/output_{args.outfile}/energies.npy',energies)
	np.savetxt(f'../test_data/output_{args.outfile}/samples_inf.txt',(args.L,args.Nc,args.ns,args.nb,args.nw,args.b))

	t4 = time.time()

	print(f'Completed. \n{n_samples} samples collected \
		\nTotal time = {t4-t0}\nWarmup time = {t2-t1}\nSampling time = {t3-t2}')

	# Close connection to samples file
	del samples

	# Show plot of energies if desired
	if args.sp:	
		import matplotlib.pyplot as plt
		plt.plot(energies)
		plt.title(f'Ising energies of samples')
		plt.show()

if __name__ == '__main__':
	main()
