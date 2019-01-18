import numpy as np

class mf_potts():
	def __init__(data):
		self.data = data
	def train(self):
		# Select training data
		n_training_samples = len(self.data.ix_train)
		training_samples = data[self.data.ix_train,:,:]
		# One-point frequencies
		f1s_train = np.mean(training_samples,axis=0)
		# Two-point frequencies
		f2s_train = np.mean([np.outer(x,x) for x in training_samples])
		# Covariance
		cov = f2s_train - np.outer(f1s_train,f1s_train)
		# Mean-field couplings
		self.J = np.linalg.inv(cov).reshape(self.data.n_spins,self.data.n_states,self.data.n_spins,self.data.n_states)
		# Set self-couplings to zero
		for i in range(self.n_spins):
			self.J[i,:,i,:] = 0
		# mean-field fields
		self.h = np.log(f1s / np.clip((1 - np.sum(f1s,axis=1)),0.001,None).reshape(self.n_spins,1)) - np.tensordot(self.J,f1s,axes=2)
	def get_fields():
	    return self.h.reshape(self.n_spins*self.n_states)
    def get_couplings():
        return self.J.reshape(self.n_spins*self.n_states,self.n_spins*self.n_states)
    def energies(seqs):
        # Field-related energies
        fE = np.tensordot(rand_samples,self.h,axes=2)
        # Coupling-related energies
        cE = np.einsum('ijk,jklm,ilm->i',rand_samples,self.J,rand_samples)/2
        return fE + cE



