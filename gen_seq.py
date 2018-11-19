import numpy as np
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of flag for magnetisation and correlation files')
parser.add_argument('-m','--model-type', type=str, action='store', dest='model_type',default='ising_mf', help='Type of model to be trained')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()

if args.model_type == 'ising_mf':
	
	# Load magnetisations and correlations from file
	#h_mf = np.genfromtxt('../test_data/output_{}/fields_{}.out'.format(args.infile,args.infile))
	#J_mf = np.genfromtxt('../test_data/output_{}/couplings_{}.out'.format(args.infile,args.infile))

	h_mf = np.zeros((40,3)).flatten()
	J_mf = np.zeros((120,120))
	J_mf[0,119] = 100
	


	def ProposeChanges(N,L):
		which_pos = np.random.randint(L,size=N)
		old_bases = samples[range(N),which_pos,:]
		new_bases = 1*(np.random.randint(4,size=(N,1))==np.arange(1,4))
		return (which_pos,old_bases,new_bases)

	def IsingEnergyChanges(samples,changes,h,J):
		(which_pos,old_bases,new_bases) = changes
		field_energy_change = np.sum(h_mf[which_pos,:] * old_bases, axis=1) - np.sum(h_mf[which_pos,:] * new_bases, axis=1)
		coupling_energy_change = np.sum(np.sum(J_mf[which_pos,:,:,:] * old_bases[...,np.newaxis,np.newaxis], axis=1) * samples, axis=(1,2)) \
			- np.sum(np.sum(J_mf[which_pos,:,:,:] * new_bases[...,np.newaxis,np.newaxis], axis=1) * samples, axis=(1,2))
		return field_energy_change + coupling_energy_change

	def MetropolisChanges(samples,deltas,changes,beta,N):
		(which_pos,old_bases,new_bases) = changes
		accepted = np.exp(-deltas*beta)>np.random.rand(N)
		accepted_changes = np.where(np.repeat(accepted.reshape(N,1),3,axis=1),new_bases,old_bases)
		if not np.sum(deltas)==0:
			print(np.sum(deltas))
		new_samples = np.copy(samples)
		new_samples[range(N),which_pos,:] = accepted_changes
		return new_samples

	def IsingEnergy(samples,h,J):
		"""Finds the Ising energy for each sample sequence"""
		field_energy = -np.tensordot(h,samples,axes=((0,1),(1,2)))
		coupling_energy = np.diag(-np.tensordot(samples,np.tensordot(samples,J,axes=((1,2),(2,3))),axes=((1,2),(1,2))))
		return field_energy + coupling_energy


	# Number of samples 
	N = 1
	L = len(h_mf)//3
	# Reshape model parameters
	h_mf = h_mf.reshape(L,3)
	J_mf = J_mf.reshape(L,3,L,3)
	# Generate random sequences to start with
	samples = 1*(np.random.randint(4,size=(N,L,1))==np.arange(1,4))
	print('Generated samples... \nOptimising samples...')
	# Minimise energy of samples by flipping one spin at a time
	n_steps = 2000
	beta = 1
	energies = np.zeros((N,n_steps))


	for j in range(n_steps):
		# propose changes

		changes = ProposeChanges(N,L)
		deltas = IsingEnergyChanges(samples,changes,h_mf,J_mf)
		samples = MetropolisChanges(samples,deltas,changes,beta,N)

		# Monitor progress
		if (j%(n_steps//10))==0:
			print('{0}/{1} steps complete'.format(j,n_steps))
		energies[:,j] = IsingEnergy(samples,h_mf,J_mf)

	# Save final samples and energies across timesteps
	np.savetxt('../test_data/samples.out',samples.reshape(N,3*L))
	np.savetxt('../test_data/energies.out',energies)

	def decode(line):
		if np.nonzero(line)[0].size:
			return np.nonzero(line)[0] + 1
		else:
			return 0 
	samples = np.apply_along_axis(decode,2,samples)
	samples_text = np.apply_along_axis(''.join,1,np.array(['A','C','G','T'])[samples])
	sample_seqs = (SeqRecord(Seq(x[0]),id=str(n)) for (n,x) in enumerate(samples_text))
	SeqIO.write(sample_seqs,'../test_data/samples_gen_lT.fasta','fasta')

	plt.plot(energies.T)
	plt.show()






