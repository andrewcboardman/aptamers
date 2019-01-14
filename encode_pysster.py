from pysster.Data import Data

def encode(filenames,structures_provided=False):
	if not structures_provided:
		data = Data(filenames,'ACGT')
	else:
		data = Data(filenames,('ACGT','FTHIMS'))
	n_spins = data.data[0].shape[0]
	n_states = data.data[0].shape[1]
	n_samples = len(data.data)
	params = (n_spins,n_states,n_samples)
	return (data,params)