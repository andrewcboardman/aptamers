import numpy as np
# Load magnetisations and correlations from file

mag = np.genfromtxt('../test_data/magnetisations.out').reshape(40,3)
corr = np.genfromtxt('../test_data/correlations.out')

# Use mean-field approximation to get model parameters
# Pseudocounts to prevent intense negative correlations from dominating

J_mf = -np.linalg.inv(corr).reshape(40,3,40,3)
h_mf = np.arctanh(mag) - np.tensordot(J_mf,mag,axes=((2,3),(0,1)))
L = len(h_mf)

# Ising 'energy' function
def IsingEnergy(seq,h,J):
	field_energy = -np.tensordot(h,seq)
	coupling_energy = -np.tensordot(seq,np.tensordot(J,seq,axes=((2,3),(0,1))))
	return field_energy + coupling_energy

# Neighbour-generating function for Gibbs sampling
def Neighbour(seq):
	new_seq  = np.array(seq)
	which_pos = np.random.randint(L)
	which_base = 1*(np.random.randint(4)==np.arange(1,4))
	new_seq[which_pos,:] = which_base
	return new_seq

# Metropolis sampler
def Metropolis(E1,E2):
	if E1 > E2:
		return 1
	else:
		#print('Energy difference of {} therefore probability of {}'.format(E1-E2,np.exp(E1-E2)))
		return 0#np.exp(E1-E2)

# Use Monte Carlo simulated annealing to draw samples with a low Ising energy

# Generate random sequences 
init_seqs = 1*(np.random.randint(3,size=(1000,L,1))==np.arange(1,4))

# Find low energy sequences in random pool to start with
#init_energies = np.zeros(1000)
#for i in range(1000):
#	init_energies[i] = IsingEnergy(init_seqs[i,...],h_mf,J_mf)
#cutoff = np.sort(init_energies)[100]
#init_seqs = init_seqs[init_energies<cutoff,...]



# Lower temperature while drawing and recording samples using Gibbs sampling
n_samples = 2000
samples = np.zeros((100,n_samples,L,3))
energies = np.zeros((100,n_samples))
for i in range(100):
	seq = init_seqs[i,...]
	for j in range(n_samples):
		beta = j/(n_samples*5)
		energy = IsingEnergy(seq,h_mf*beta,J_mf*beta)
		new_seq = Neighbour(seq)
		new_energy = IsingEnergy(new_seq,h_mf*beta,J_mf*beta)
		if Metropolis(energy,new_energy) > np.random.rand():
			print(energy-new_energy)
			seq = new_seq
			energy = new_energy
			
		samples[i,j,...] = seq
		energies[i,j] = energy/beta

import matplotlib.pyplot as plt
(fig,(ax1,ax2)) = plt.subplots(1,2)
ax1.plot(energies.T)
ax2.imshow(np.sum(samples[:,1999,...],axis=0))

plt.show()

new_mag = np.sum(samples[:,1999,...],axis=0)