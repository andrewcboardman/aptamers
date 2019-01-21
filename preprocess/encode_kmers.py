import numpy as np
import argparse
from Bio import SeqIO
import itertools
from Bio import SeqRecord

def kmer_content(record,k,kmer_hash):
	string = str(record.seq)
	kmer_content = np.unique([kmer_hash[string[i:i+k]] for i in range(len(string)-k)],return_counts=True)
	return kmer_content

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Flag for original model')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
	parser.add_argument('-k', type=int,action='store')
	args = parser.parse_args()

	# Generate kmer hash
	bases = ('A','C','G','T')
	kmers = [''.join(x) for x in itertools.product(bases,repeat=args.k)]
	kmer_hash = dict(zip(kmers,range(4**args.k)))

	# Read FASTA-formatted samples
	samples = SeqIO.parse(f'../test_data/output_{args.infile}/{args.outfile}/samples.fasta','fasta')

	# Calculate kmer counts
	kmer_contents = (kmer_content(record,args.k,kmer_hash) for record in samples)

	# Write to file
	with open(f'../test_data/output_{args.infile}/{args.outfile}/samples_{args.k}mers.txt','w') as file:
		for record in kmer_contents:
			file.write(' '.join(record[0].astype(str)) + '\n' + ' '.join(record[1].astype(str)) + '\n')

if __name__ == '__main__':
	main()