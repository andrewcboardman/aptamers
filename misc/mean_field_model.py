import numpy as np

# Load magnetisations and correlations from file

mag = np.genfromtxt('../test_data/magnetisations.out')
corr = np.genfromtxt('../test_data/correlations.out')

# Calculate parameters of model

J_mf = np.linalg.inv(corr)
h_mf = np.arctanh(mag) - np.dot(J_mf,mag)

# Save to file

np.savetxt('../test_data/couplings.out',J_mf)
np.savetxt('../test_data/fields.out',h_mf)
