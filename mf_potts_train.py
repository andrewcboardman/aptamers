import argparse
from encode_pysster import encode
from mf_potts import mf_potts
from histo_energies import plot

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Input filename')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='infile',	help='Output filename')
args = parser.parse_args()

data, params = encode(args.infile)
model = mf_potts(data,params)
mf_potts.data.train_val_test_split(0.9,0.1)

print('Imported data. Infering model...')
J_mf,h_mf = model.train()

print('Inference complete, saving coefficients...')
np.savetxt(args.outfile,J_mf)

test_energies,random_energies = model.test()
plot(test_energies,random_energies)
