import numpy as np

# Load magnetisations and correlations from file

mag = np.genfromtxt('magnetisations.out')
corr = np.genfromtxt('correlations.out')

# Calculate parameters of model

J_mf = np.linalg.inverse(corr)
h_mf = np.arctanh(mag) - np.dot(J_mf,mag)

# Save to file

np.savetxt(J_mf,'couplings.out')
np.savetxt(h_mf,'fields.out')
