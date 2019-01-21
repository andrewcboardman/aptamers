import argparse
from pysster.Data import Data
from pysster.Model import Model
from pysster import utils
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i1', '--infile1', type=str, action='store', dest='infile1',	help='Input filename 1')
parser.add_argument('-i2', '--infile2', type=str, action='store', dest='infile2',	help='Input filename')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',	help='Output folder')
args = parser.parse_args()

output_folder = args.outfile
if not os.path.exists(output_folder):
	os.mkdir(output_folder)

data = Data([args.infile1,args.infile2],('ACGT','FHIMST'))
data.train_val_test_split(0.8,0.1,1)

model_params = {'conv_num':2, 'kernel_len':20}
model = Model(model_params,data)
model.print_summary()

print('Imported data. Infering model...')
model.train(data)

predictions = model.predict(data,'test')
labels = data.get_labels('test')
utils.plot_roc(labels, predictions, output_folder+"roc.png")
utils.plot_prec_recall(labels, predictions, output_folder+"prec.png")
print(utils.get_performance_report(labels, predictions))
activations = model.get_max_activations(data, "test")
logos = model.visualize_all_kernels(activations, data, output_folder)
model.plot_clustering(activations, output_folder+"clustering.png")

utils.save_data(data, output_folder+"data.pkl")
utils.save_model(model, output_folder+"model.pkl")