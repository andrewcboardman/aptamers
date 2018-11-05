import numpy as np
from Bio import SeqIO
import pickle
def encode(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# convert into a data vector
	boole = 1*(chars[:,np.newaxis] == ['A','C','G']).flatten()
	# Join vector to string
	return np.array2string(boole)

def SeqEncoder(seqfile):
	# Import reads from file
	with open('output.txt','w') as fid:
		for seqline in SeqIO.parse(seqfile,'fastq'):
			# encode each read as a binary data vector
			codedline = encode(str(seqline))
			fid.write(codedline+'\n')
