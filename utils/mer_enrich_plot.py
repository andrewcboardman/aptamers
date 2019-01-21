import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
k = 9
plus_recs = list(SeqIO.parse(f'plus_mers/mers_{k}.fasta','fasta'))
plus_seqs = [str(x.seq) for x in plus_recs]
plus_counts = [int(x.id) for x in plus_recs]
minus_recs = list(SeqIO.parse(f'minus_mers/mers_{k}.fasta','fasta'))
minus_seqs = [str(x.seq) for x in minus_recs]
minus_counts = [int(x.id) for x in minus_recs]
all_seqs = np.unique(plus_seqs + minus_seqs)

all_plus_counts = np.zeros_like(all_seqs,dtype=int)
all_minus_counts = np.zeros_like(all_seqs,dtype=int)
for i,seq in enumerate(all_seqs):
	if seq in plus_seqs:
		all_plus_counts[i] = plus_counts[plus_seqs.index(seq)]
	if seq in minus_seqs:
		all_minus_counts[i] = minus_counts[minus_seqs.index(seq)]

plt.scatter(all_plus_counts,all_minus_counts,color='red',alpha=0.5)
plt.plot([0,max(all_plus_counts)],[0,max(all_minus_counts)],color='black')
plt.xlim([0,max(all_plus_counts)])
plt.ylim([0,max(all_minus_counts)])
plt.xlabel(f'Number of occurences of {k}mer in the binding sequences')
plt.ylabel(f'Number of occurences of {k}mer in the control sequences')
plt.show()