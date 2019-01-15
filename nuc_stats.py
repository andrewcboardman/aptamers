from pysster.Data import Data
import matplotlib.pyplot as plt
data = Data(['minus.fasta'],'ACGT')
profile = sum((data.astype('int32') for data in data.data))/len(data.data)
plt.switch_backend('TkAgg')
plt.bar(range(40),profile[:,0])
plt.bar(range(40),profile[:,1],bottom=profile[:,0])
plt.bar(range(40),profile[:,2],bottom=profile[:,0]+profile[:,1])
plt.bar(range(40),profile[:,3],bottom=profile[:,0]+profile[:,1]+profile[:,2])
plt.legend(('A','C','G','T')
	)
plt.show()