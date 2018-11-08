import numpy as np

J = np.genfromtxt('couplings_rand.out')
h = np.genfromtxt('fields_rand.out')

def Energy(seq):
	field_energy = np.dot(h,seq)
	coupling_energy = np.dot(seq,np.dot(J,seq))
	return field_energy + coupling_energy


# Generate a random initial bit sequence
seq = np.random.choice((-1,1),size=120)

# generate a variety of bit sequences using Gibbs sampling and average
steps = 1000
seq_energy = Energy(seq)
wu_energies = np.zeros(steps)
for i in range(steps):
	pos = np.random.randint(0,120)
	newseq = seq
	if seq[pos] == 1:
		newseq[pos] = 0
	else:
		newseq[pos] = 1
	new_energy = Energy(newseq)
	if (new_energy < seq_energy) or (np.random.rand() < np.exp((seq_energy-new_energy)/0.1)):
		seq = newseq
		seq_energy = Energy(seq)
	wu_energies[i] = seq_energy


store_seqs = np.zeros((steps,120))
store_energies = np.zeros(steps)
count_stored = 0

for i in range(steps):
	pos = np.random.randint(0,120)
	newseq = seq
	if seq[pos] == 1:
		newseq[pos] = 0
	else:
		newseq[pos] = 1
	if (Energy(newseq) < seq_energy) or (np.random.rand() < np.exp((seq_energy-new_energy)/0.01)):
		seq = newseq
		seq_energy = Energy(seq)
	store_seqs[count_stored,:] = seq
	store_energies[count_stored] = seq_energy
	count_stored += 1

store_seqs = store_seqs[:count_stored]
store_energies = store_energies[:count_stored]

print(np.tensordot(store_seqs,store_seqs,axes=(0,0)))
import seaborn
import matplotlib.pyplot as plt
seaborn.heatmap(np.tensordot(store_seqs,store_seqs,axes=(0,0)))
plt.show()
# compare base enrichment at each site
