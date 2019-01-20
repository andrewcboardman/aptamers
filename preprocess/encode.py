import numpy as np
import argparse
import os
from Bio import SeqIO
import itertools
import load_samples

def seq2code(seq):
	seq_arr = np.array(list(seq))[:,np.newaxis]
	code_arr = np.stack(3*[np.zeros_like(seq_arr,dtype=int)])
	array[seqarr==np.array(('C','G','T'))] = 1
	return array

# def seq2array_indels(seq,n_spins):
# 	array = np.zeros((n_spins,4),dtype=int)
# 	seqarr = np.array(list(seq))[:,np.newaxis]
# 	array[seqarr==np.array(('A','C','G','T'))] = 1
# 	return array


# def encode_kmers(record,k,kmer_hash):
# 	string = str(record.seq)
# 	kmer_content = np.unique([kmer_hash[string[i:i+k]] for i in range(len(string)-k)],return_counts=True)
# 	kmer_spins = np.zeros(4**k,dtype=int)
# 	kmer_spins[kmer_content[0]] = 1
# 	return kmer_spins

# def seq2array_ortho(seq,n_spins):
# 	array = np.zeros((n_spins,3),dtype=int)
# 	for i in range(n_spins):
# 		if seq[i] == 'A':
# 			array[i,:] = (-3,0,0)
# 		elif seq[i] == 'C':
# 			array[i,:] = (1,-2,0)
# 		elif seq[i] == 'G':
# 			array[i,:] = (1,1,-1)
# 		elif seq[i] == 'T':
# 			array[i,:] = (1,1,1)
# 	return array



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='FASTA-formatted samples')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='name of encoded text file')
	parser.add_argument('-m', '--mode', type=str, action='store', dest='mode', default='bases',help='Method of encoding')
	parser.add_argument('-k', type=int,action='store')
	args = parser.parse_args()

	# Read FASTA-formatted samples
	seqs = SeqIO.parse(args.path + args.infile,'fasta')


	# load sample metadata
	n_samples,n_spins,b = load_samples.load_metadata(args.path+'samples_inf.txt')

	if args.mode == 'bases':
		strings = (str(seq.seq) for seq in seqs)
		coded_strings = (seq2code(string) for string in strings)
		with open(args.path + args.outfile,'w') as file:
			for record in coded_strings:
				file.write(' '.join(record.astype(str).flatten()))

	# elif args.mode == 'bases_align':
	# 	strings = (str(seq.seq) for seq in seqs)
	# 	coded_strings = (seq2array_indels(string,n_spins) for string in strings)
	# 	with open(args.path + args.outfile,'w') as file:
	# 		for record in coded_strings:
	# 			file.write(' '.join(record[:,0].astype(str)) + '\n' + \
	# 				' '.join(record[:,1].astype(str)) + '\n' + \
	# 				' '.join(record[:,2].astype(str)) + '\n' + \
	# 				' '.join(record[:,3].astype(str)) + '\n')

	# elif args.mode == 'kmers':
	# 	# Generate kmer hash
	# 	bases = ('A','C','G','T')
	# 	kmers = [''.join(x) for x in itertools.product(bases,repeat=args.k)]
	# 	kmer_hash = dict(zip(kmers,range(4**args.k)))

	# 	# Calculate kmer counts
	# 	kmer_contents = (encode_kmers(record,args.k,kmer_hash) for record in seqs)

	# 	with open(args.path + args.outfile,'w') as file:
	# 		for record in kmer_contents:
	# 			file.write(' '.join(record.astype(str)) + '\n')

	# elif args.mode == 'bases_ortho':
	# 	strings = (str(seq.seq) for seq in seqs)
	# 	coded_strings = (seq2array_ortho(string,n_spins) for string in strings)
	# 	with open(args.path + args.outfile,'w') as file:
	# 		for record in coded_strings:
	# 			file.write(' '.join(record[:,0].astype(str)) + '\n' + \
	# 				' '.join(record[:,1].astype(str)) + '\n' + \
	# 				' '.join(record[:,2].astype(str)) + '\n')

if __name__ == '__main__':
	main()
