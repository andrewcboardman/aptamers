from pysster.Data import Data
import numpy as np
import matplotlib.pyplot as plt
data = Data(['plus.fasta'],'ACGT')
covar = sum((np.outer(data.astype('int32'),data.astype('int32')) for data in data.data)).reshape(40,4,40,4)/len(data.data)
plt.switch_backend('TkAgg')
fig, axs = plt.subplots(3,3)
for i in range(3):
	for j in range(3):
		covar_tmp = covar[:,i,:,j]
		covar_tmp[range(40),range(40)] = np.mean(covar_tmp)
		axs[i,j].imshow(covar_tmp)
plt.title('Plots of C, G and T covariance ')
plt.show()