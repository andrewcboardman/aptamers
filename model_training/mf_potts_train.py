import argparse
from Data import Data
from mf_potts import mf_potts
from histo_energies import plot
from scikit_learn import metrics

parser = argparse.ArgumentParser()
parser.add_argument('-i1', '--infile1', type=str, action='store', dest='infile1',	help='Input filename 1')
parser.add_argument('-i2', '--infile2', type=str, action='store', dest='infile2',	help='Input filename 2')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='outfile',	help='Output filename')
args = parser.parse_args()

print('Importing data...')
data_plus = Data('args.infile1',40,3)
data_plus.split(0.8,0.0,1)
data_minus = Data('args.infile1',40,3)
data_minus.split(0.8,0.0,1)
model_plus = mf_potts(data_plus)
model_minus = mf_potts(data_minus)

print('Imported data. Inferring model...')

model_plus.train()
J_plus = model_plus.get_couplings()
h_plus = model_plus.get_fields()

model_minus.train()
J_minus = model_minus.get_couplings()
h_minus = model_minus.get_fields()

print('Inference complete, testing...')

plus_energies_train_plus = model_plus.energies(data_plus[data_plus.ix_train])
minus_energies_train_plus = model_minus.energies(data_plus[data_plus.ix_train])
P_train_plus = np.exp(- plus_energies_train_plus/minus_energies_train_plus)

plus_energies_train_minus = model_plus.energies(data_minus[data_minus.ix_train])
minus_energies_train_minus = model_minus.energies(data_minus[data_minus.ix_train])
P_train_minus = np.exp(- plus_energies_train_minus/minus_energies_train_minus)

train_true = np.concatenate((np.ones_like(P_train_plus),np.zeros_like(P_train_minus))
train_prob = np.concatenate(P_train_plus,P_train_minus)
train_auc = metrics.roc_auc_score(train_true,train_prob)
print(f'Training AUC {train_auc}')

plus_energies_test_plus = model_plus.energies(data_plus[data_plus.ix_test])
minus_energies_test_plus = model_minus.energies(data_plus[data_plus.ix_test])
P_test_plus = np.exp(- plus_energies_test_plus/minus_energies_test_plus)

plus_energies_test_minus = model_plus.energies(data_minus[data_minus.ix_test])
minus_energies_test_minus = model_minus.energies(data_minus[data_minus.ix_test])
P_test_minus = np.exp(- plus_energies_test_minus/minus_energies_test_minus)

test_true = np.concatenate((np.ones_like(P_test_plus),np.zeros_like(P_test_minus))
test_prob = np.concatenate(P_test_plus,P_test_minus)
test_auc = metrics.roc_auc_score(test_true,test_prob)
print(f'Test AUC {test_auc}')
