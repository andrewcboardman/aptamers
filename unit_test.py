import numpy as np

J_init = np.random.rand(120,120)
h_init = np.zeros(120)


def ProposeChanges(N,L,samples):
		which_pos = np.random.randint(L,size=N)
		old_bases = samples[range(N),which_pos,:]
		new_bases = 1*(np.random.randint(4,size=(N,1))==np.arange(1,4))
		base_changes = new_bases - old_bases
		return (which_pos,base_changes)

def IsingEnergyChanges(samples,changes,h,J):
	(which_pos,base_changes) = changes
	field_energy_changes = -np.sum(h[which_pos,:] * base_changes, axis=1) 
	# Coupling energy change = product of base changes with current samples across J (which is symmetric in leading diagonal)
	coupling_energy_changes = -2*np.sum(np.sum(J[which_pos,:,:,:] * base_changes[...,np.newaxis,np.newaxis], axis=1) * samples, axis=(1,2)) 
	return field_energy_changes + coupling_energy_changes

def MetropolisChanges(samples,deltas,changes,beta,N):
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


# Number of samples 
N = 100
L = len(h_init)//3
# Reshape model parameters
h_init = h_init.reshape(L,3)
J_init = J_init.reshape(L,3,L,3)
# Generate random sequences to start with
init_samples = 1*(np.random.randint(4,size=(N,L,1))==np.arange(1,4))

print('Generated samples... \nOptimising samples...')
# Minimise energy of samples by flipping one spin at a time
n_steps = 200
beta = 1
energies = np.zeros((N,n_steps))

samples = np.zeros((n_steps+1,N,L,3))
samples[0,:,:,:] = init_samples

for j in range(n_steps):
	# propose changes

	changes = ProposeChanges(N,L,samples[j,:,:,:])
	deltas = IsingEnergyChanges(samples[j,:,:,:],changes,h_init,J_init)
	samples[j+1,:,:,:] = MetropolisChanges(samples[j,:,:,:],deltas,changes,beta,N)

	# Monitor progress
	if (j%(n_steps//10))==0:
		print('{0}/{1} steps complete'.format(j,n_steps))
	energies[:,j] = IsingEnergy(samples[j,:,:,:],h_init,J_init)


def decode(line):
	if np.nonzero(line)[0].size:
		return np.nonzero(line)[0] + 1
	else:
		return 0 

samples = np.apply_along_axis(decode,3,samples)
samples_text = np.apply_along_axis(''.join,2,np.array(['A','C','G','T'])[samples])
samples_text = samples_text.flatten(samples_text.shape[0]*samples_text.shape[1])

def parse_infile(infile):
	if infile.endswith('fasta'):
		return SeqIO.parse(infile,'fasta')
		print('Reading {}...'.format(infile))
	elif infile.endswith('fastq'):
		return SeqIO.parse(infile,'fastq')
		print('Reading {}...'.format(infile))
	else:
		print('Warning: unrecognised input file format')
def averages(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# One-hot encoding to find magnetisations
	mag1s = np.where(chars.reshape(L,1) == ['C','G','T'],1,0).flatten()
	# two-site magnetisations from outer product
	mag2s = np.outer(mag1s,mag1s)
	return (mag1s,mag2s)

mag1s_av = np.zeros(3*L)
mag2s_av = np.zeros((3*L,3*L))
print('Counting frequencies...')
mags = (averages(seq) for seq in samples_text)
for mag in mags:
	mag1s_av += mag[0]/N
	mag2s_av += mag[1]/N

# correlations = average two-site magnetisations - outer product of one-site magnetisations
corr = mag2s_av - np.outer(mag1s_av,mag1s_av)
mag = mag1s_av
print('Calculated frequencies.')

J_mf = -np.linalg.inv(corr)
for i in range(120):
	J_mf[i:(i+3),i:(i+3)] = 0
h_mf = np.arctanh(mag) - np.dot(J_mf,mag)
import matplotlib.pyplot as plt
(fig,axs) = plt.subplots(1,2)

im2 = axs[1].imshow(J_mf)
fig.colorbar(im1)
fig.colorbar(im2)
plt.show()