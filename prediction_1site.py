import numpy as np

# Load magnetisations from file

mag = np.genfromtxt('../test_data/magnetisations.out')

# Process magnetisation vector to generate a profile matrix

L = int(len(mag)/4)
mag_rs = np.reshape(mag,(4,L),order='F')
p_blank = np.abs(1 - np.sum(mag_rs,axis=0))
profile = np.vstack((mag_rs,p_blank))

# Generate sequences matching this profile

seq = np.empty(L,dtype=str)
for i in range(L):
	seq[i] = np.random.choice(['A','C','G','T','-'],p=profile[:,i])
print(''.join(seq))