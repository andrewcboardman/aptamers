import numpy as np
# In the Curie-Weiss model, all the couplings are the same
L = 40 
J_cw = np.zeros((L,3,L,3))
h_cw = np.zeros((L,3))
h_cw[:,0] = 1


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
N = 100
n_steps = 4000
# Generate random sequences to start with
init_samples = 1*(np.random.randint(4,size=(N,L,1))==np.arange(1,4))
samples = np.zeros((n_steps,N,L,3))

energies = np.zeros((N,n_steps))
print('Generated samples... \n')

print('Optimising samples...')
# Minimise energy of samples by flipping one spin at a time

beta = 1

# Warmup markov chain
for j in range(n_steps//2):
	# Calculate energy
	energies[:,j] = IsingEnergy(samples[j,...],h_cw,J_cw)
	# Generate new samples and calculate energy
	new_samples = Neighbours(samples[j,...])
	new_energies = IsingEnergy(new_samples,h_cw,J_cw)
	# Change samples if necessary
	metro = Metropolis(energies[:,j],new_energies,beta)
	samples[j,...] = np.where(metro,new_samples,samples[j,...])
	# Monitor progress
	if j==n_steps//2:
		print('{0}/{1} steps complete'.format(j,n_steps))

# Now run markov chain and store values

for j in range(n_steps//2,n_steps):
	# Calculate energy
	energies[:,j] = IsingEnergy(init_samples,h_cw,J_cw)
	# Generate new samples and calculate energy
	new_samples = Neighbours(init_samples)
	new_energies = IsingEnergy(new_samples,h_cw,J_cw)
	# Change samples if necessary
	metro = Metropolis(energies[:,j],new_energies,beta)
	init_samples = np.where(metro,new_samples,samples[j,...])


# Save final samples and energies across timesteps
np.savetxt('../test_data/samples.out',samples.reshape(n_steps*N,3*L))
np.savetxt('../test_data/energies.out',energies)


import matplotlib.pyplot as plt
(fig1,(ax1,ax2)) = plt.subplots(1,2)
ax1.plot(energies.T)
ax2.imshow(np.sum(samples[n_steps//2:,...],axis=(0,1)))

plt.show()