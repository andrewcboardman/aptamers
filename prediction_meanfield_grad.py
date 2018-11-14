
import numpy as np
# Load magnetisations and correlations from file

mag = np.genfromtxt('../test_data/magnetisations.out').reshape(40,3)
corr = np.genfromtxt('../test_data/correlations.out')
# Use mean-field approximation to get model parameters
# Pseudocounts to prevent intense negative correlations from dominating
# Reshape to 
J_mf = -np.linalg.inv(corr).reshape(40,3,40,3)
h_mf = np.arctanh(mag) - np.tensordot(J_mf,mag,axes=((2,3),(0,1)))
L = len(h_mf)

def IsingEnergy(samples,h,J):
	"""Finds the Ising energy for each sample sequence"""
	field_energy = -np.tensordot(h,samples,axes=((0,1),(1,2)))
	coupling_energy = np.diag(-np.tensordot(samples,np.tensordot(samples,J,axes=((1,2),(2,3))),axes=((1,2),(1,2))))
	return field_energy + coupling_energy

def Neighbours(samples):
	"""Generates a new batch of sequences
	 which differ from the current ones by one base each"""
	new_samples  = np.array(samples)
	# Select a base in each sequence to modify
	which_pos = np.random.randint(L,size=N)
	which_base = 1*(np.random.randint(4,size=(N,1))==np.arange(1,4))
	new_samples[:,which_pos,:] = which_base
	return new_samples

def Metropolis(energy,new_energy,beta):
	delta = (-energy+new_energy)
	metro = np.exp(-delta*beta)>np.random.rand(*delta.shape)
	return metro[:,np.newaxis,np.newaxis]


# Number of samples 
N = 1000
# Generate random sequences to start with
samples = 1*(np.random.randint(3,size=(N,L,1))==np.arange(1,4))

# Minimise energy of samples by flipping one spin at a time
n_steps = 500
beta = 0
energies = np.zeros((n_steps,N))
for j in range(n_steps):
	# Calculate energy
	energies[j,:] = IsingEnergy(samples,h_mf,J_mf)
	# Generate new samples and calculate energy
	new_samples = Neighbours(samples)
	new_energies = IsingEnergy(new_samples,h_mf,J_mf)
	# Change samples if necessary
	metro = Metropolis(energies[j,:],new_energies,beta)
	samples = np.where(metro,samples,new_samples)

import matplotlib.pyplot as plt
(fig,(ax1,ax2)) = plt.subplots(1,2)
ax1.plot(energies.T)
ax2.imshow(np.sum(samples[:,1999,...],axis=0))

plt.show()

new_mag = np.sum(samples[:,1999,...],axis=0)