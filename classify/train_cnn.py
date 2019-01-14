from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.layers import Embedding
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
import numpy as np

seq_length = 40
x_train = np.random.randint(2,size=(1000, 40, 3))
y_train = np.random.randint(2, size=(1000, 1))
x_test = np.random.randint(2,size=(100, 40, 3))
y_test = np.random.randint(2, size=(100, 1))

model = Sequential()
model.add(Conv1D(40, 3, activation='relu', input_shape=(seq_length, 3)))
model.add(Conv1D(40, 3, activation='relu'))
model.add(MaxPooling1D(3))
model.add(Conv1D(80, 3, activation='relu'))
model.add(Conv1D(80, 3, activation='relu'))
model.add(GlobalAveragePooling1D())
model.add(Dropout(0.5))
model.add(Dense(1, activation='sigmoid'))

model.compile(loss='binary_crossentropy',
              optimizer='rmsprop',
              metrics=['accuracy'])

model.fit(x_train, y_train, batch_size=50, epochs=10)
score = model.evaluate(x_test, y_test, batch_size=50)
