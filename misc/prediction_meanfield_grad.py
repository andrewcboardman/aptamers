
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of flag for magnetisation and correlation files')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()

# Load magnetisations and correlations from file

#mag = np.genfromtxt('../test_data/magnetisations_{}.out'.format(args.infile))+0.01
#corr = np.genfromtxt('../test_data/correlations_{}.out'.format(args.infile)) 
#mag = mag.reshape(mag.shape[0]//3,3)

# Use mean-field approximation to get model parameters
# Pseudocounts to prevent intense negative correlations from dominating
# Reshape to 
#J_mf = -np.linalg.inv(corr).reshape(corr.shape[0]//3,3,corr.shape[1]//3,3)
#h_mf = np.arctanh(mag) - np.tensordot(J_mf,mag,axes=((2,3),(0,1)))
#L = len(h_mf)


h_mf = np.zeros((40,3))
J_mf = np.zeros((120,120))
J_mf[0,119] = 100
J_mf= J_mf.reshape(40,3,40,3)
L = 40
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
samples = 1*(np.random.randint(4,size=(N,L,1))==np.arange(1,4))
print('Generated samples... \nOptimising samples...')
# Minimise energy of samples by flipping one spin at a time
n_steps = 2000
beta = 0.01
energies = np.zeros((N,n_steps))
for j in range(n_steps):
	# Calculate energy
	energies[:,j] = IsingEnergy(samples,h_mf,J_mf)
	# Generate new samples and calculate energy
	new_samples = Neighbours(samples)
	new_energies = IsingEnergy(new_samples,h_mf,J_mf)
	# Change samples if necessary
	metro = Metropolis(energies[:,j],new_energies,beta)
	samples = np.where(metro,new_samples,samples)
	# Monitor progress
	if j==n_steps//2:
		print('{0}/{1} steps complete'.format(j,n_steps))

# Save final samples and energies across timesteps
np.savetxt('../test_data/samples.out',samples.reshape(N,3*L))
np.savetxt('../test_data/energies.out',energies)

# plo

import matplotlib.pyplot as plt
(fig,(ax1,ax2)) = plt.subplots(1,2)
ax1.plot(energies.T)
ax2.imshow(np.sum(samples,axis=0))

plt.show()
