import numpy as np
from Bio import SeqIO

input_file = 'test.fastq'
score_cutoff = 15


# Count number of reads
count = 0
for rec in SeqIO.parse(input_file, "fastq"):
	count += 1
print('{} total reads'.format(count))

# Filter for any bases which have quality lower than a cutoff
good_reads = (rec for rec in \
	SeqIO.parse(input_file, "fastq") \
	if min(rec.letter_annotations["phred_quality"]) >= score_cutoff)
count = SeqIO.write(good_reads, "test_filter.fastq", "fastq")
print("Saved {} reads".format(count))

# 


def SeqEncoder(seqfile):
	# Import reads from file
	with open('output.txt','w') as fid:
		for seqline in SeqIO.parse(seqfile,'fastq'):

			# encode each read as a binary data vector
			codedline = encode(str(seqline))
			fid.write(codedline+'\n')



def encode(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# convert into a data vector
	boole = 1*(chars[:,np.newaxis] == ['A','C','G']).flatten()
	# Join vector to string
	return np.array2string(boole)


