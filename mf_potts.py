import numpy as np

class mf_potts():
	def __init__(data,params):
		self.data = data
		self.n_spins = params[0]
		self.n_states = params[1]
		self.n_samples = params[2]
	def train(self):
		# Select training data
		n_training_samples = len(data._get_idx('train'))
		training_samples = [data.data[i] for i in data._get_idx('train')]
		# One-point frequencies
		f1s_train = sum(training_samples)/n_training_samples
		# Two-point frequencies
		f2s_train = sum((np.outer(x,x) for x in training_samples))/n_training_samples
		# Covariance
		cov = f2s_train - np.outer(f1s_train,f1s_train)
		# Mean-field couplings
		self.J = np.linalg.inv(cov).reshape(self.n_spins,self.n_states,self.n_spins,self.n_states)
		# Set self-couplings to zero
		for i in range(self.n_spins):
			self.J[i,:,i,:] = 0
		# mean-field fields
		self.h = np.log(f1s / np.clip((1 - np.sum(f1s,axis=1)),0.001,None).reshape(self.n_spins,1)) - np.tensordot(self.J,f1s,axes=2)
		return self.J.reshape(self.n_spins*self.n_states,self.n_spins*self.n_states),self.h
	def test(self):
		# Select validation data
		val_samples = np.stack([data.data[i] for i in data._get_idx('val')])
		# Generate set of random sequences of same size
		random_samples = np.zeros_like(val_samples)
		random_samples[:,:,np.random.randint(self.n_states,size=self.n_spins)] = 1
		# Field energies for validation samples
		val_energies = np.tensordot(val_samples,self.h,axes=2)
		# Coupling energies for validation samples
		val_energies -= np.einsum('ijk,jklm,ilm->i',val_samples,self.J,val_samples)/2
		# Field energies for random samples
		rand_energies = np.tensordot(rand_samples,self.h,axes=2)
		# Coupling energies for random samples
		rand_energies -= np.einsum('ijk,jklm,ilm->i',rand_samples,self.J,rand_samples)/2
		return (val_energies,rand_energies)



