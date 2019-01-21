from Bio import SeqIO
from Bio import SeqRecord
import os
import numpy as np
import argparse
def random_base(n):
	return ''.join(np.array(('A','C','G','T'))[np.random.randint(4,size=abs(n))])
def shift(seqrecord,n_shift):
	if n_shift == 0:
		return seqrecord
	elif n_shift > 0:
		return random_base(n_shift) + seqrecord[:-n_shift]
	else:
		return seqrecord[(-n_shift):] + random_base((n_shift))

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Input file in FASTA format')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Output file in FASTA format')
parser.add_argument('-s', '--shift_length', type=int, action='store', dest='shift_length', help='Number of bases to shift samples by')
args = parser.parse_args()

input_seqs = SeqIO.parse(args.infile,'fasta')
output_seqs = (shift(seq,args.shift_length) for seq in input_seqs)
with open(args.outfile, "w") as output_handle:
    SeqIO.write(output_seqs, output_handle, "fasta")

