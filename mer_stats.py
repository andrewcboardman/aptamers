import numpy as np
import matplotlib.pyplot as plt
unique = np.zeros(30)
distinct = np.zeros(30)
total = np.zeros(30)
maxcount = np.zeros(30)
for k in range(3,30):
	with open(f'mers_{k+1}.stats','r') as file:
		rawtext = file.readlines()
	unique[k] = int(rawtext[0].split()[-1])
	distinct[k] = int(rawtext[1].split()[-1])
	total[k] = int(rawtext[2].split()[-1])
	maxcount[k] = int(rawtext[3].split()[-1])

fig, axs = plt.subplots(2,2)
axs[0,0].plot(unique)
axs[0,0].set_title('No. k-mers found only once')
axs[1,0].plot(distinct)
axs[1,0].set_title('No. distinct k-mers found')
axs[0,1].plot(total)
axs[0,1].set_title('Total no. k-mers found')
axs[1,1].plot(maxcount)
axs[1,1].set_title('No. of counts of most frequent k-mer')
plt.show()