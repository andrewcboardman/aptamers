from Bio import SeqIO
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
args = parser.parse_args()

# Read FASTA-formatted samples
seqs = SeqIO.parse(f'../test_data/output_{args.infile}/{args.outfile}/samples.fasta','fasta')	
# Create numpy array (stored on hard drive to prevent RAM overloading) 
samples = np.memmap(f'../test_data/output_{args.infile}/{args.outfile}/samples.npy',mode='w+',dtype=int,shape=(Nc*ns//nw,L,3))
# Populate numpy array with encoded samples
encode.encode(seqs,samples,L)

f1s = np.zeros((L,3))
f2s = np.zeros((L*3,L*3))
print('Counting frequencies...')
fs = ((seq,np.outer(seq,seq)) for seq in samples)
for f in fs:
	f1s += f[0]/N
	f2s += f[1]/N

f1s = np.concatenate((f1s,(1 - np.sum(f1s,axis=1)).reshape(L,1)),axis=1)
f2s = f2s.reshape(L,3,L,3)
f2s = np.concatenate((f2s,(1 - np.sum(f2s,axis=1)).reshape(L,1,L,3)),axis=1)
f2s = np.concatenate((f2s,(1 - np.sum(f2s,axis=3)).reshape(L,3,L,1)),axis=3)

fig, (ax1,ax2) = plt.subplots(1,2)
ax1.plot(f1s)
ax2.plot(f2s)
plt.show()