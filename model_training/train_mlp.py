from keras import Sequential

def main():

	data = # <- load standardised data
	model = Sequential()# <-Some list of layers
	model.compile() # <- training specification
	model.train()

	# Test model and save outputs