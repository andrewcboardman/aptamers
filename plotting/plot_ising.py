import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Folder containing original couplings')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
parser.add_argument('-s', '--slice_length',type=int, action='store', dest='slice_length',default=0, help='length of slice from sequence used in training')
args = parser.parse_args()

J_init = np.load(f'../test_data/output_{args.infile}/couplings_init.npy')
h_init = np.load(f'../test_data/output_{args.infile}/fields_init.npy')
if args.slice_length == 0:
	J_mf = np.load(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf.npy')
	h_mf = np.load(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf.npy')

	fig, (ax1,ax2) = plt.subplots(1,2)

	ax1.scatter(h_init.flatten(),h_mf.flatten(),color='red',alpha=0.5)
	ax1.legend()
	ax1.set_xlabel(r'h$_{init}$')
	ax1.set_ylabel(r'h$_{MF}$')


	ax2.scatter(J_init.flatten(),J_mf.flatten(),color='red',alpha=0.5)
	ax2.legend()
	ax2.set_xlabel(r'J$_{init}$')
	ax2.set_ylabel(r'J$_{MF}$')
	print(np.corrcoef(J_init.flatten(),J_mf.flatten()))

	plt.suptitle('Actual vs. inferred values for mean-field inference of a Gaussian matrix padded with zeros')
	plt.show()
else:
	J_mf_slice = np.load(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf_slice_{args.slice_length}.npy')
	h_mf_slice = np.load(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf_slice_{args.slice_length}.npy')

	L = J_init.shape[0]
	gap = (L - args.slice_length) //2
	J_mf = np.zeros((L,3,L,3))
	J_mf[gap:(L-gap),:,gap:(L-gap),:] = J_mf_slice
	h_mf = np.zeros((L,3))
	h_mf[gap:(L-gap),:] = h_mf_slice
	fig, (ax1,ax2) = plt.subplots(1,2)

	ax1.scatter(h_init.flatten(),h_mf.flatten(),color='red',alpha=0.5)
	ax1.set_xlabel(r'h$_{init}$')
	ax1.set_ylabel(r'h$_{MF}$')


	ax2.scatter(J_init.flatten(),J_mf.flatten(),color='red',alpha=0.5)
	ax2.set_xlabel(r'J$_{init}$')
	ax2.set_ylabel(r'J$_{MF}$')
	print(np.corrcoef(J_init.flatten(),J_mf.flatten()))

	plt.suptitle('Actual vs. inferred values for mean-field inference of a Gaussian matrix padded with zeros')
	plt.show()

#fit2 = np.polyfit(J_init.flatten(),J_mf.flatten(),1)
#ax2.plot(J_init.flatten(),fit2[0]*J_init.flatten() + fit2[1])


#J_init_std = np.std(J_init)
#J_mf_std = np.std(J_mf)
#corr = np.corrcoef(J_init.flatten(),J_mf.flatten())[0,1]
#ax2.text(0,0,r'$\sigma$(J$_{MF}$)'+f'= {J_mf_std:.2f}')
#ax2.text(0,0.5,r'$\sigma$(J$_{init}$)'+f'= {J_init_std:.2f}')
#ax2.text(0,0,r'Correlation '+f'= {h_mf_std:.2f}')
#im = ax3.imshow(np.triu(J_init.reshape(60,60))+np.tril(J_mf.reshape(60,60)))
#ax3.set_xlabel(r'J$_{init}$')
#ax3.set_ylabel(r'J$_{MF}$')
#plt.colorbar(im)




