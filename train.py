
import numpy as np
import argparse
from Bio import SeqIO
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',default='test.fastq',	help='Name of flag for magnetisation and correlation files')
parser.add_argument('-m','--model-type', type=str, action='store', dest='model_type',default='ising_mf', help='Type of model to be trained')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',default='', help='Flag for output files')
args = parser.parse_args()


def parse_infile(infile):
	if args.infile.endswith('fasta'):
		return SeqIO.parse(args.infile,'fasta')
		print('Reading {}...'.format(args.infile))
	elif args.infile.endswith('fastq'):
		return SeqIO.parse(args.infile,'fastq')
		print('Reading {}...'.format(args.infile))
	else:
		print('Warning: unrecognised input file format')

# Check pre-processed sequence data in FASTA or FASTQ format
# Reads must be all the same length and not contain ambiguous base calls
seqs = parse_infile(args.infile)
N = 0
L = 0
for seq in seqs:
	N += 1
	L = max(L,len(seq))
print('Found {} reads, of length {}'.format(N,L))

# Set up folder for output
if not os.path.isdir('../test_data/output_{}'.format(args.outfile)):
	os.mkdir('../test_data/output_{}'.format(args.outfile))

# Choose which model to train

if args.model_type == 'ising_mf':
	# Ising model with mean-field approximation
	# Parse sequences and extract <s_i> and <s_i*s_j> (One- and two-site frequencies)
	def averages(seq):
		# Split up sequence into chars
		chars = np.array(list(seq),dtype=str)
		# One-hot encoding to find magnetisations
		mag1s = np.where(chars.reshape(L,1) == ['C','G','T'],1,0).flatten()
		# two-site magnetisations from outer product
		mag2s = np.outer(mag1s,mag1s)
		return (mag1s,mag2s)

	mag1s_av = np.zeros(3*L)
	mag2s_av = np.zeros((3*L,3*L))
	print('Counting frequencies...')
	seqs = parse_infile(args.infile)
	mags = (averages(seq) for seq in seqs)
	for mag in mags:
		mag1s_av += mag[0]/N
		mag2s_av += mag[1]/N

	# correlations = average two-site magnetisations - outer product of one-site magnetisations
	corr = mag2s_av - np.outer(mag1s_av,mag1s_av)
	mag = mag1s_av
	print('Calculated frequencies.')

	# Save magnetisations and correlations
	np.savetxt('../test_data/output_{}/magnetisations_{}.out'.format(args.outfile,args.outfile),mag1s_av)
	np.savetxt('../test_data/output_{}/correlations_{}.out'.format(args.outfile,args.outfile),corr)

	# Add pseudocounts to regularise
	#mag = (3/4)* mag + 1/4
	#corr = (3/4)*corr + 1/4

	# Use mean-field approximation to find fields and couplings
	print('Calculating fields and couplings')
	import matplotlib.pyplot as plt
	im = plt.imshow(corr)
	plt.colorbar(im)
	plt.show()
	J_mf = -np.linalg.inv(corr)
	for i in range(120):
		J_mf[i:(i+3),i:(i+3)] = 0
	h_mf = np.arctanh(mag) - np.dot(J_mf,mag)

	# Save fields and couplings
	np.savetxt('../test_data/output_{}/fields_{}.out'.format(args.outfile,args.outfile),h_mf)
	np.savetxt('../test_data/output_{}/couplings_{}.out'.format(args.outfile,args.outfile),J_mf)

elif args.model_type == 'ising_ace':
	print("This hasn't been implemented yet lol")


