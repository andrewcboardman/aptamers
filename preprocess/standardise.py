import numpy as np

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='FASTA-formatted samples')
	parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='name of encoded text file')
	args = parser.parse_args()

	# Read one-hot encoded file
    data = np.genfromtxt(args.infile)

    # Calculate means and standard deviations
    mean = np.mean(data, axis=0)
    stdev = np.std(data,axis=0)
    standardised_data = (data - mean) / stdev

    # Write to file
    np.savetxt(args.outfile,standardised_data)

if __name__ == '__main__':
	main()
