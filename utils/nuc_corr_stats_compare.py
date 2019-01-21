from pysster.Data import Data
import matplotlib.pyplot as plt
import numpy as np
data_minus = Data(['minus.fasta'],'ACGT')
data_plus = Data(['plus.fasta'],'ACGT')
f2s_plus = sum((np.outer(data.astype('int32'),data.astype('int32')) for data in data_plus.data)).reshape(40,4,40,4)/len(data_plus.data)
f2s_minus = sum((np.outer(data.astype('int32'),data.astype('int32')) for data in data_minus.data)).reshape(40,4,40,4)/len(data_minus.data)
plt.switch_backend('TkAgg')
bases = ['A','C','G','T']
for i in range(4):
	for j in range(4):
		plt.scatter(f2s_plus[:,i,:,j].flatten(),f2s_minus[:,i,:,j].flatten(),label=f'{bases[i]}-{bases[j]}')
plt.plot([0,0.15],[0,0.15])
plt.title('Comparison of two-site frequencies across the binding and control datasets')
plt.xlabel('Frequency of base in binding sequence')
plt.ylabel('Frequency of base in control sequence')
plt.xlim([0,0.15])
plt.ylim([0,0.15])
plt.legend()
plt.show()