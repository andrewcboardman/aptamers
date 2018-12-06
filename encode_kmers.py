import numpy as np
import argparse
from Bio import SeqIO
import itertools
from Bio import SeqRecord

def annotate(record,k,kmer_hash):
	string = str(record.seq)
	kmer_content = [str(kmer_hash[string[i:i+k]]) for i in range(len(string)-k)]
	return SeqRecord.SeqRecord(record.seq,id=record.id,description = ' '.join(kmer_content))		

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

	# Add k-mer counts as descriptions
	samples_annotated = (annotate(record,args.k,kmer_hash) for record in samples)

	SeqIO.write(samples_annotated,f'../test_data/output_{args.infile}/{args.outfile}/samples_{args.k}mers.fasta','fasta')

if __name__ == '__main__':
	main()