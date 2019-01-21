import argparse
from encode_pysster import encode
from mf_potts import mf_potts
from histo_energies import plot
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Input filename (binding sequences)')
parser.add_argument('-c', '--control', type=str, action='store', dest='control',	help='Input filename (control sequences)')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',	help='Output filename')
args = parser.parse_args()

data, params = encode([args.infile])
data_c, params_c = encode([args.control])
model = mf_potts(data,params)
model.data.train_val_test_split(0.9,0.1)
data_c.train_val_test_split(0.9,0.1)
print('Imported data. Inferring model coefficients...')
J_mf,h_mf = model.train()

print('Inference complete, saving coefficients...')
np.savetxt(args.outfile,np.concatenate((J_mf,h_mf),axis=0))

test_energies,control_energies,random_energies = model.test(data_c)

plot(test_energies,control_energies)
