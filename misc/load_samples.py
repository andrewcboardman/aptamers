import numpy as np
def load_samples(filename,n_samples,n_spins,n_states):
	raw_samples = np.genfromtxt(filename,dtype=int)
	samples = raw_samples.reshape(n_samples,n_states,n_spins)
	return samples.transpose(0,2,1)

def load_metadata(filename):
	metadata = np.genfromtxt(filename,dtype=float)
	n_spins, Nc, ns, nb, nw = metadata[:-1].astype('int32')
	b = metadata[-1]
	n_samples = Nc*ns//nw
	return (n_samples,n_spins,b)
	
