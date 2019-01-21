import numpy as np
import time
class mf_potts():
	def __init__(self,data):
		self.data = data
	def train(self):
		t0 = time.time()
		# Select training data
		n_training_samples = len(self.data.ix_train)
		training_samples = np.squeeze(self.data.data[self.data.ix_train,:,:])
		# One-point frequencies
		f1s = np.mean(training_samples,axis=0)
		# Two-point frequencies
		f2s = sum((np.outer(x,x) for x in training_samples))/n_training_samples
		# Covariance
		cov = f2s - np.outer(f1s,f1s)
		# Mean-field couplings
		self.J = -np.linalg.inv(cov).reshape(self.data.n_spins,self.data.n_states,self.data.n_spins,self.data.n_states)
		# Set self-couplings to zero
		for i in range(self.data.n_spins):
			self.J[i,:,i,:] = 0
		# mean-field fields
		self.h = np.log(f1s / np.clip((1 - np.sum(f1s,axis=1)),0.001,None).reshape(self.data.n_spins,1)) - np.tensordot(self.J,f1s,axes=2)
	def get_fields(self):
		return self.h.reshape(self.data.n_spins*self.data.n_states)
	def get_couplings(self):
		return self.J.reshape(self.data.n_spins*self.data.n_states,self.data.n_spins*self.data.n_states)
	def energies(self,samples):
		# Field-related energies
		fE = -np.tensordot(samples,self.h,axes=2)
		# Coupling-related energies
		cE = -np.einsum('ijk,jklm,ilm->i',samples,self.J,samples)/2
		return fE + cE



