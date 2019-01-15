from sklearn import svm
import numpy as np
filenames = ['../real_data/data/minus.fasta','../real_data/data/plus.fasta']
from pysster.Data import Data
data = Data(filenames,'ACGT')
data.train_val_test_split(0.8,0.0,1)
X_train = np.stack([data.data[i].flatten() for i in data._get_idx('train')])
print(X_train.shape)
y_train = data.get_labels('train')[:,1]

clf = svm.SVC(gamma='scale')
clf.fit(X_train, y_train)