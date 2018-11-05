import numpy as np
from Bio import SeqIO
import pandas
input_file = 'test.fastq'
score_cutoff = 15


# Count number of reads

count = 0
maxlen = 0
for rec in SeqIO.parse(input_file, "fastq"):
	count += 1
	if len(rec) > maxlen:
		maxlen = len(rec)
print('{} total reads'.format(count))

# Filter for any reads which have minimum base quality lower than a cutoff
### Filter for ambiguous calls?

good_reads = (rec for rec in \
	SeqIO.parse(input_file, "fastq") \
	if min(rec.letter_annotations["phred_quality"]) >= score_cutoff)
count = SeqIO.write(good_reads, "test_filter.fastq", "fastq")
print("Saved {} reads".format(count))

# Encode reads in 1-of-3 format with adenine as baseline

def encode(seq):
	# Split up sequence into chars
	chars = np.array(list(seq),dtype=str)
	# convert into a data vector
	boole = 1*(chars[:,np.newaxis] == ['C','G','T']).flatten()
	# Join vector to string
	return boole

# write encoded reads to file

output_data = np.zeros(count,)
for seqline in good_reads:
		# encode each read as a binary data vector
		codedline = encode(str(seqline))
		fid.write(codedline+'\n')






