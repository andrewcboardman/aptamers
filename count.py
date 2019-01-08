import numpy as np


def load_samples(filename,n_samples,n_spins,n_states):
	raw_samples = np.genfromtxt(filename,dtype=int)
	samples = raw_samples.reshape(n_samples,n_states,n_spins)
	return samples.transpose(0,2,1)

def count(samples,n_samples,n_spins,n_states):
	f1s = np.zeros((n_spins,n_states))
	f2s = np.zeros((n_spins*n_states,n_spins*n_states))
	print('Counting frequencies...')
	fs = ((seq,np.outer(seq,seq)) for seq in samples)
	for f in fs:
		f1s += f[0]/n_samples
		f2s += f[1]/n_samples
	return (f1s,f2s)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m','--model', type=str, action='store', dest='model',default='mean_field', help='Type of model to be trained')
	parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')

	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='File containing encoded samples')

	parser.add_argument('-N', '--n_samples', type=str, action='store', help='Number of samples')
	parser.add_argument('-q', '--n_states', type=str, action='store', help='Number of states per spin')
	parser.add_argument('-L', '--n_spins', type=str, action='store', help='Number of spins')

	parser.add_argument('-c', '--correct', action='store_true', dest='correct', help='Apply correction to fields')

	parser.add_argument('-s', '--slice_length', type=int, action='store', dest='slice_length', help='Length of slice to take of sample')
	args = parser.parse_args()



	samples = load_samples(args.path + args.infile,args.n_samples,args.n_spins,args.n_states)

	f1s, f2s = count(samples,args.n_samples,args.n_spins,args.n_states)

	np.savetxt(args.path+'f1s_'+args.infile,f1s)
	np.savetxt(args.path+'f2s_'+args.infile,f2s)

