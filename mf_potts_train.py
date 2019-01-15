import argparse
from encode_pysster import encode
from mf_potts import mf_potts
from histo_energies import plot
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Input filename')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',	help='Output filename')
args = parser.parse_args()

data, params = encode([args.infile])
model = mf_potts(data,params)
model.data.train_val_test_split(0.9,0.1)

print('Imported data. Inferring model coefficients...')
J_mf,h_mf = model.train()

print('Inference complete, saving coefficients...')
np.savetxt(args.outfile,np.concatenate((J_mf,h_mf),axis=0))

test_energies,random_energies = model.test()

print(test_energies.shape,random_energies.shape)
print(test_energies)
print(random_energies)
plot(test_energies,random_energies)
