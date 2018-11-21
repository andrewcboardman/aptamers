from sklearn import svm
import numpy as np

x_train = np.array(((0,0),(1,0),(0,1),(1,1),(4,4)))
x_test = np.array(((1,0),(1,1)))
y_train = np.array((0,1,1,0,0))
y_test = np.array((1,0))

clf = svm.SVC(gamma='scale',kernel='rbf')
clf.fit(x_train,y_train)
print('predictions are',clf.predict(x_test))
print('support vectors are',clf.support_vectors_)

import matplotlib.pyplot as plt
plt.scatter(x_train[:,0],x_train[:,1],c=y_train)