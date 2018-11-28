import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Folder containing original couplings')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile', help='Flag for output files')
args = parser.parse_args()

J_init = np.load(f'../test_data/output_{args.infile}/couplings_init.npy')
h_init = np.load(f'../test_data/output_{args.infile}/fields_init.npy')

J_mf = np.load(f'../test_data/output_{args.infile}/{args.outfile}/couplings_mf.npy')
h_mf = np.load(f'../test_data/output_{args.infile}/{args.outfile}/fields_mf.npy')

fig, (ax1,ax2) = plt.subplots(1,2)


ax1.scatter(h_init,h_mf,color='red',alpha=0.5)
ax1.set_xlabel(r'h$_{init}$')
ax1.set_ylabel(r'h$_{MF}$')
#ax1.text(0,0,r'$\sigma$(h$_{MF}$)'+f'= {h_mf_std:.2f}')

ax2.scatter(J_init,J_mf,color='red',alpha=0.05)
ax2.set_xlabel(r'J$_{init}$')
ax2.set_ylabel(r'J$_{MF}$')
#fit2 = np.polyfit(J_init.flatten(),J_mf.flatten(),1)
#ax2.plot(J_init.flatten(),fit2[0]*J_init.flatten() + fit2[1])


#J_init_std = np.std(J_init)
#J_mf_std = np.std(J_mf)
#corr = np.corrcoef(J_init.flatten(),J_mf.flatten())[0,1]
#ax2.text(0,0,r'$\sigma$(J$_{MF}$)'+f'= {J_mf_std:.2f}')
#ax2.text(0,0.5,r'$\sigma$(J$_{init}$)'+f'= {J_init_std:.2f}')
#ax2.text(0,0,r'Correlation '+f'= {h_mf_std:.2f}')
plt.suptitle('Actual vs. inferred values for mean-field inference of a frameshifted Gaussian matrix')
#im = ax3.imshow(np.triu(J_init.reshape(60,60))+np.tril(J_mf.reshape(60,60)))
#ax3.set_xlabel(r'J$_{init}$')
#ax3.set_ylabel(r'J$_{MF}$')
#plt.colorbar(im)

plt.show()


