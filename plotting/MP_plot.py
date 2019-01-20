import matplotlib.pyplot as plt
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-p','--path',type=str,action='store', dest='path',help='path to working folder')
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='File containing samples')
parser.add_argument('-N', '--n_samples', type=int, action='store', help='Number of samples')
parser.add_argument('-q', '--n_states', type=int, action='store', help='Number of states per spin')
parser.add_argument('-L', '--n_spins', type=int, action='store', help='Number of spins')
args = parser.parse_args()

# load 1 and 2-site frequencies
f1s = np.genfromtxt(args.path+'f1s_' + args.infile).flatten()
f2s = np.genfromtxt(args.path+'f2s_' + args.infile)#.reshape(args.n_spins,args.n_states,args.n_spins,args.n_states)

# Calculate covariance matrix
corr = f2s - np.outer(f1s,f1s)#.reshape(args.n_spins,args.n_states,args.n_spins,args.n_states)
eigvals = np.zeros((args.n_states*args.n_spins,args.n_states*args.n_spins))

# Calculate eigenvalue spectrum
#for i in range(args.n_states):
#	for j in range(args.n_states):
#		eigvals[i,j,:] = np.linalg.eigvals(corr[:,i,:,j])
eigvals = np.linalg.eigvals(corr)
# Largest eigenvalue (used as a scale for graph)
rad = np.max(eigvals)

# Fit Marchenko-Pastur distribution
g = args.n_spins/args.n_samples
x = np.linspace(0.0001,rad,10000)
y = np.sqrt(np.abs((1 + np.sqrt(g))**2 - x)*np.abs(x - (1 - np.sqrt(g))**2))/(2*np.pi*x)

# Plot eigenvalue histogram against MP distribution
plt.hist(eigvals.flatten(),bins=50,density=True)
plt.plot(x,y)
plt.xlabel('Magnitude of eigenvalue')
plt.ylabel('Probability')
plt.show()
