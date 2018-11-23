import numpy as np
#############################################################################################

# Generating sequences from an ising model
J_init = np.random.rand(120,120)
J_init = 0.5*(J_init+J_init.T)
h_init = np.zeros(120)

# Number of samples 
N = 1000
L = len(h_init)//3
# Reshape model parameters
h_init = h_init.reshape(L,3)
J_init = J_init.reshape(L,3,L,3)
# Generate random sequences to start with
init_samples = 1*(np.random.randint(4,size=(N,L,1))==np.arange(1,4))

print('Generated samples... \nOptimising samples...')
# Minimise energy of samples by flipping one spin at a time
n_steps = 500
beta = 1
energies = np.zeros((N,10))

samples = np.zeros((n_steps+1,N,L,3))
samples[0,:,:,:] = init_samples



def ProposeChanges(N,L,samples):
	# Randomly generate proposed changes to the samples
	which_pos = np.random.randint(L,size=N)
	old_bases = samples[range(N),which_pos,:]
	new_bases = 1*(np.random.randint(4,size=(N,1))==np.arange(1,4))
	base_changes = new_bases - old_bases
	return (which_pos,base_changes)

def IsingEnergyChanges(samples,changes,h,J):
	""" Calculate the changes in energy associated with each of the proposed changes"""
	(which_pos,base_changes) = changes
	field_energy_changes = -np.sum(h[which_pos,:] * base_changes, axis=1) 
	# Coupling energy change = product of base changes with current samples across J (which is symmetric in leading diagonal)
	coupling_energy_changes = -2*np.sum(np.sum(J[which_pos,:,:,:] * base_changes[...,np.newaxis,np.newaxis], axis=1) * samples, axis=(1,2)) 
	return field_energy_changes + coupling_energy_changes

def MetropolisChanges(samples,deltas,changes,beta,N):
	"""Accept or reject changes based on their energy change"""
	(which_pos,base_changes) = changes
	accepted = np.exp(-deltas*beta)>np.random.rand(N)
	accepted_no = np.arange(N)[accepted]
	accepted_pos = which_pos[accepted]
	accepted_changes = base_changes[accepted]
	new_samples = np.copy(samples)
	new_samples[accepted_no,accepted_pos,:] += accepted_changes
	return new_samples

def IsingEnergy(samples,h,J):
	"""Finds the Ising energy for each sample sequence"""
	field_energy = -np.tensordot(h,samples,axes=((0,1),(1,2)))
	coupling_energy = np.diag(-np.tensordot(samples,np.tensordot(samples,J,axes=((1,2),(2,3))),axes=((1,2),(1,2))))
	return field_energy + coupling_energy

# Run Monte Carlo to generate samples from Ising model
for j in range(n_steps):
	changes = ProposeChanges(N,L,samples[j,:,:,:])
	deltas = IsingEnergyChanges(samples[j,:,:,:],changes,h_init,J_init)
	samples[j+1,:,:,:] = MetropolisChanges(samples[j,:,:,:],deltas,changes,beta,N)
	# Monitor progress
	if (j%(n_steps//10))==0:
		print('{0}/{1} steps complete'.format(j,n_steps))
		energies[:,j//(n_steps//10)] = IsingEnergy(samples[j,:,:,:],h_init,J_init)

# select last set of samples which should have relaxed to Boltzmann distribution
samples = samples[-1,...].reshape(N,40,3)


# Calculate one and two-point frequencies
f1s = np.zeros(3*L)
f2s = np.zeros((3*L,3*L))
print('Counting frequencies...')
for sample in samples:
	f1s += sample.flatten()
	f2s += np.outer(sample,sample)
f1s = (f1s + 0/4)/(1*N)
f2s = (f2s + 0/4**2)/(1*N)

# Calculate covariances
cov = f2s - np.outer(f1s,f1s)
print('Calculated frequencies.')

print('Calculating couplings...')
J_mf = -np.linalg.inv(cov) 
for i in range(120):
	J_mf[i:(i+3),i:(i+3)] = 0
h_mf = np.arctanh(f1s) - np.dot(J_mf,f1s)


import matplotlib.pyplot as plt
fig,axs = plt.subplots(1,3)
axs[0].plot(energies.T)
im1 = axs[1].imshow(np.tril(J_mf) + np.triu(J_init.reshape(3*L,3*L)),interpolation='none')
im2 = axs[2].plot(np.stack((h_init.reshape(120),h_mf)).T)
plt.colorbar(im1)
plt.show()

#fig,ax = plt.subplots(1,1)
#ax.hist(energies[:,-1])
#minE = np.min(energies[:,-1])
#maxE = np.max(energies[:,-1])
#ax.plot(np.linspace(minE,maxE,100),np.exp(-np.linspace(0,maxE-minE,100)))
#plt.show()