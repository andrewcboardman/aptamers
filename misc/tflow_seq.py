import numpy as np
from tensorflow import keras

train_x_partial = x[1000:,...]
train_x_val = x[1000:1500,...]
t


model = keras.Sequential([
    keras.layers.Flatten(input_shape=(40, 3)),
    keras.layers.Dense(128, activation=tf.nn.relu),
    keras.layers.Dense(10, activation=tf.nn.softmax)])


model.compile(optimizer=tf.train.AdamOptimizer(), 
              loss='binary_crossentropy',
              metrics=['accuracy'])

model.fit(train_data, train_labels, epochs=5)

results = model.evaluate(test_data, test_labels)



