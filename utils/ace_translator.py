import numpy as np
f1s = np.genfromtxt('../test_data/output_test40/1m/f1s_samples_bases.txt')
f2s = np.genfromtxt('../test_data/output_test40/1m/f2s_samples_bases.txt').reshape(40,3,40,3)
with open('samples_count.p','w') as file:
	for i in range(40):
		file.write(' '.join(f1s[i,:].astype(str))+'\n')
	for j in range(40):
		for k in range(j+1,40):
			file.write(' '.join(f2s[j,:,k,:].astype(str).flatten())+'\n')
