import argparse
from encode_pysster import encode
from pysster.Model import Model
from pysster import utils

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Input filename')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='infile',	help='Output filename')
args = parser.parse_args()

data, params = encode(args.infile)
data.train_val_test_split(0.8,0.1,0.1)

model_params = {'conv_num':1, 'kernel_len':10}
model = Model(data,model_params)

print('Imported data. Infering model...')
model.train(data)

predictions = model.predict(data,'test')
utils.plot_roc(labels, predictions, output_folder+"roc.png")
utils.plot_prec_recall(labels, predictions, output_folder+"prec.png")
print(utils.get_performance_report(labels, predictions))
activations = model.get_max_activations(data, "test")
logos = model.visualize_all_kernels(activations, data, output_folder)

utils.save_as_meme([logo[0] for logo in logos], output_folder+"motifs_seq.meme")
utils.save_as_meme([logo[1] for logo in logos], output_folder+"motifs_struct.meme")


model.plot_clustering(activations, output_folder+"clustering.png")
utils.save_data(data, output_folder+"data.pkl")
utils.save_model(model, output_folder+"model.pkl")
