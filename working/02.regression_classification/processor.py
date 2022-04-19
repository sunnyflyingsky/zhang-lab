"""
# Author: Yuhan Fei
# Created Time :  15 May 2020
# File Name: loader.py
# Description: Data loader
"""

import numpy as np
import pandas as pd

import pickle
import os
import pathlib

this_dir = str(pathlib.Path(__file__).parent.absolute())

MAX_ATOM = 400
MAX_BOND = MAX_ATOM * 2


def roc_curve(y_pred, y_label, figure_file, method_name):
    '''
		y_pred is a list of length n.  (0,1)
		y_label is a list of same length. 0/1
		https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html#sphx-glr-auto-examples-model-selection-plot-roc-py  
	'''
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc
    from sklearn.metrics import roc_auc_score
    y_label = np.array(y_label)
    y_pred = np.array(y_pred)
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    fpr[0], tpr[0], _ = roc_curve(y_label, y_pred)
    roc_auc[0] = auc(fpr[0], tpr[0])
    lw = 2
    plt.plot(fpr[0], tpr[0],
             lw=lw, label=method_name + ' (area = %0.2f)' % roc_auc[0])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    fontsize = 14
    plt.xlabel('False Positive Rate', fontsize=fontsize)
    plt.ylabel('True Positive Rate', fontsize=fontsize)
    plt.title('Receiver Operating Characteristic Curve')
    plt.legend(loc="lower right")
    plt.savefig(figure_file)
    return

def prauc_curve(y_pred, y_label, figure_file, method_name):
    '''
		y_pred is a list of length n.  (0,1)
		y_label is a list of same length. 0/1
		reference: 
			https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
	'''
    import matplotlib.pyplot as plt
    from sklearn.metrics import precision_recall_curve, average_precision_score
    from sklearn.metrics import f1_score
    from sklearn.metrics import auc
    lr_precision, lr_recall, _ = precision_recall_curve(y_label, y_pred)
    #	plt.plot([0,1], [no_skill, no_skill], linestyle='--')
    plt.plot(lr_recall, lr_precision, lw=2,
             label=method_name + ' (area = %0.2f)' % average_precision_score(y_label, y_pred))
    fontsize = 14
    plt.xlabel('Recall', fontsize=fontsize)
    plt.ylabel('Precision', fontsize=fontsize)
    plt.title('Precision Recall Curve')
    plt.legend()
    plt.savefig(figure_file)
    return

def length_func(list_or_tensor):
    if type(list_or_tensor) == list:
        return len(list_or_tensor)
    return list_or_tensor.shape[0]

def index_select_ND(source, dim, index):
    index_size = index.size()
    suffix_dim = source.size()[1:]
    final_size = index_size + suffix_dim
    target = source.index_select(dim, index.view(-1))
    return target.view(final_size)

"""
import networkx as nx
def atom_features(atom):
	return \
		np.array(one_of_k_encoding_unk(atom.GetSymbol(),
									   ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'As',
										'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti',
										'Zn', 'H', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt',
										'Hg', 'Pb', 'Unknown']) +  # 44
				 one_of_k_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) +  # 获取原子连接数（受H是否隐藏影响）11
				 one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) +  # 与该原子连接的氢原子个数 11
				 one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) +  # 获取原子隐式化合价 11
				 [atom.GetIsAromatic()])  # 1

def one_of_k_encoding(x, allowable_set):
	if x not in allowable_set:
		raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
	return list(map(lambda s: x == s, allowable_set))

def one_of_k_encoding_unk(x, allowable_set):
	if x not in allowable_set:
		x = allowable_set[-1]
	return list(map(lambda s: x == s, allowable_set))

from rdkit import Chem
def smiles2graph(smile):
	mol = Chem.MolFromSmiles(smile)

	c_size = mol.GetNumAtoms()

	features = []
	for atom in mol.GetAtoms():
		feature = atom_features(atom)
		features.append(feature / sum(feature))

	edges = []
	for bond in mol.GetBonds():
		edges.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
	g = nx.Graph(edges).to_directed()
	edge_index = []
	for e1, e2 in g.edges:
		edge_index.append([e1, e2])

	return c_size, features, edge_index
"""

"""
from dgllife.utils import smiles_to_bigraph
def smiles2dgl_canonical(s, node_featurizer, edge_featurizer):
	try:
		g = smiles_to_bigraph(smiles=s,
                  node_featurizer=node_featurizer,
                  edge_featurizer=edge_featurizer)
	except:
		print('dgl canonical fingerprint not working for smiles: ' + s + ' removed...')
		g = 'REMOVED'
	return g
"""

##############################################
## 针对drug, target(protein/RNA) sequence, target strcture, target knots进行编码 
##############################################
def encode_drug(df_data, drug_encoding, column_name='SMILES', save_column_name='drug_encoding'):
    print('encoding drug...')
    print('unique drugs: ' + str(len(df_data[column_name].unique())))

    if drug_encoding == 'CNN':
        unique = pd.Series(df_data[column_name].unique()).apply(trans_drug)
        unique_dict = dict(zip(df_data[column_name].unique(), unique))
        df_data[save_column_name] = [unique_dict[i] for i in df_data[column_name]]
    # the embedding is large and not scalable but quick, so we move to encode in dataloader batch.
    elif drug_encoding == 'RNN':
        unique = pd.Series(df_data[column_name].unique()).apply(trans_drug)
        unique_dict = dict(zip(df_data[column_name].unique(), unique))
        df_data[save_column_name] = [unique_dict[i] for i in df_data[column_name]]
    elif drug_encoding == 'Transformer':
        unique = pd.Series(df_data[column_name].unique()).apply(drug2emb_encoder)
        unique_dict = dict(zip(df_data[column_name].unique(), unique))
        df_data[save_column_name] = [unique_dict[i] for i in df_data[column_name]]
    elif drug_encoding == 'DGL_GCN':
        df_data[save_column_name] = df_data[column_name]
    elif drug_encoding == 'DGL_GAT':
        df_data[save_column_name] = df_data[column_name]
    elif drug_encoding in ['DGL_NeuralFP']:
        df_data[save_column_name] = df_data[column_name]
    elif drug_encoding == 'DGL_AttentiveFP':
        df_data[save_column_name] = df_data[column_name]
    elif drug_encoding in ['DGL_GIN_AttrMasking', 'DGL_GIN_ContextPred']:
        df_data[save_column_name] = df_data[column_name]

    else:
        raise AttributeError("Please use the correct drug encoding available!")
    return df_data

def encode_protein(df_data, target_encoding, column_name='Target Sequence', save_column_name='target_encoding'):
    print('encoding protein...')
    print('unique target sequence: ' + str(len(df_data[column_name].unique())))

    #if target_encoding == 'CNN':
    if target_encoding == 'CNN':
        AA = pd.Series(df_data[column_name].unique()).apply(trans_protein)
        AA_dict = dict(zip(df_data[column_name].unique(), AA))
        df_data[save_column_name] = [AA_dict[i] for i in df_data[column_name]]
    # the embedding is large and not scalable but quick, so we move to encode in dataloader batch.
    elif target_encoding == 'Prism':
        AA = pd.Series(df_data[column_name].unique()).apply(trans_protein)
        AA_dict = dict(zip(df_data[column_name].unique(), AA))
        df_data[save_column_name] = [AA_dict[i] for i in df_data[column_name]]
    # the embedding is large and not scalable but quick, so we move to encode in dataloader batch.
    elif target_encoding == 'RNN':
        AA = pd.Series(df_data[column_name].unique()).apply(trans_protein)
        AA_dict = dict(zip(df_data[column_name].unique(), AA))
        df_data[save_column_name] = [AA_dict[i] for i in df_data[column_name]]

    elif target_encoding == 'Transformer':
        AA = pd.Series(df_data[column_name].unique()).apply(protein2emb_encoder)
        AA_dict = dict(zip(df_data[column_name].unique(), AA))
        df_data[save_column_name] = [AA_dict[i] for i in df_data[column_name]]
    else:
        raise AttributeError("Please use the correct protein encoding available!")
    return df_data

def encode_structure(df_data, structure_encoding, column_name='Target Structure', save_column_name='structure_encoding'):
    print('encoding structure...')
    print('unique structure sequence: ' + str(len(df_data[column_name].unique())))

    if structure_encoding == '0/1/2':
        AA = pd.Series(df_data[column_name].unique()).apply(trans_structure)
        AA_dict = dict(zip(df_data[column_name].unique(), AA))
        df_data[save_column_name] = [AA_dict[i] for i in df_data[column_name]]
    # the embedding is large and not scalable but quick, so we move to encode in dataloader batch.
    else:
        raise AttributeError("Please use the correct protein encoding available!")
    return df_data

def encode_knot(df_data, knot_encoding, column_name='Target Knot', save_column_name='knot_encoding'):
    print('encoding knot...')
    print('unique knot sequence: ' + str(len(df_data[column_name].unique())))

    if knot_encoding == '0/1/2':
        AA = pd.Series(df_data[column_name].unique()).apply(trans_knot)
        AA_dict = dict(zip(df_data[column_name].unique(), AA))
        df_data[save_column_name] = [AA_dict[i] for i in df_data[column_name]]
    # the embedding is large and not scalable but quick, so we move to encode in dataloader batch.
    else:
        raise AttributeError("Please use the correct protein encoding available!")
    return df_data

##############################################
## 组织pandas dataframe的格式，组装samples
##############################################
def split_encode(X_drug=None, X_target=None, y=None, drug_encoding=None, target_encoding=None, frac=[0.7, 0.1, 0.2],
                 random_seed=1):
    if (X_drug is not None) and (X_target is not None):
        if (X_drug is None) or (X_target is None):
            raise AttributeError("Target pair sequence should be in X_target, X_drug")

    print('Drug Target Interaction Prediction...')

    if isinstance(X_target, str):
        X_target = [X_target]
    if len(X_target) == 1:
        # one target high throughput screening setting
        X_target = np.tile(X_target, (length_func(X_drug),))

    df_data = pd.DataFrame(zip(X_drug, X_target, y))
    df_data.rename(columns={0: 'SMILES',
                            1: 'Target Sequence',
                            2: 'Label'},
                   inplace=True)
    print('in total: ' + str(len(df_data)) + ' drug-target pairs')

    df_data = encode_drug(df_data, drug_encoding)
    df_data = encode_protein(df_data, target_encoding)

    # dti split
    print('splitting dataset...')
    train_frac, val_frac, test_frac = frac
    test = df_data.sample(frac=test_frac, replace=False, random_state=random_seed)
    train_val = df_data[~df_data.index.isin(test.index)]
    val = train_val.sample(frac=val_frac / (1 - test_frac), replace=False, random_state=1)
    train = train_val[~train_val.index.isin(val.index)]

    print('Done.')
    return train.reset_index(drop=True), val.reset_index(drop=True), test.reset_index(drop=True)

def encode(X_drug=None, X_target=None, y=None, drug_encoding=None, target_encoding=None, random_seed=1):
    if (X_drug is not None) and (X_target is not None):
        if (X_drug is None) or (X_target is None):
            raise AttributeError("Target pair sequence should be in X_target, X_drug")

    print('Drug Target Interaction Prediction...')

    if isinstance(X_target, str):
        X_target = [X_target]
    if len(X_target) == 1:
        # one target high throughput screening setting
        X_target = np.tile(X_target, (length_func(X_drug),))

    df_data = pd.DataFrame(zip(X_drug, X_target, y))
    df_data.rename(columns={0: 'SMILES',
                            1: 'Target Sequence',
                            2: 'Label'},
                   inplace=True)
    print('in total: ' + str(len(df_data)) + ' drug-target pairs')

    df_data = encode_drug(df_data, drug_encoding)
    df_data = encode_protein(df_data, target_encoding)

    # dti split
    print('splitting dataset...')
    test = df_data.sample(frac=1, replace=False, random_state=random_seed)

    print('Done.')
    return test.reset_index(drop=True)

def encode_ss(X_drug=None, X_target=None, X_structure=None, y=None, drug_encoding=None, target_encoding=None, random_seed=1):
    #if (X_drug is not None) and (X_target is not None) and (X_structure is not None):
    #    if (X_drug is None) or (X_target is None) or (X_structure is not None):
    #        raise AttributeError("Target pair sequence should be in X_target, X_drug")

    print('Drug Target Interaction Prediction...')

    if isinstance(X_target, str):
        X_target = [X_target]
    if len(X_target) == 1:
        # one target high throughput screening setting
        X_target = np.tile(X_target, (length_func(X_drug),))

    df_data = pd.DataFrame(zip(X_drug, X_target, X_structure, y))
    df_data.rename(columns={0: 'SMILES',
                            1: 'Target Sequence',
                            2: 'Target Structure',
                            3: 'Label'},
                   inplace=True)
    print('in total: ' + str(len(df_data)) + ' drug-target pairs')

    df_data = encode_drug(df_data, drug_encoding)
    df_data = encode_protein(df_data, target_encoding)
    df_data = encode_structure(df_data, "0/1/2")

    # dti split
    print('splitting dataset...')
    test = df_data.sample(frac=1, replace=False, random_state=random_seed)

    print('Done.')
    return test.reset_index(drop=True)
    # 在sanmple过程中打乱了index，reset_index可以重新按照0-n的顺序建立index

def encode_ss2(X_drug=None, X_target=None, X_structure=None, X_knot=None, y=None, drug_encoding=None, target_encoding=None, random_seed=1):
    #if (X_drug is not None) and (X_target is not None) and (X_structure is not None):
    #    if (X_drug is None) or (X_target is None) or (X_structure is not None):
    #        raise AttributeError("Target pair sequence should be in X_target, X_drug")

    print('Drug Target Interaction Prediction...')

    if isinstance(X_target, str):
        X_target = [X_target]
    if len(X_target) == 1:
        # one target high throughput screening setting
        X_target = np.tile(X_target, (length_func(X_drug),))

    df_data = pd.DataFrame(zip(X_drug, X_target, X_structure, X_knot, y))
    df_data.rename(columns={0: 'SMILES',
                            1: 'Target Sequence',
                            2: 'Target Structure',
                            3: 'Target Knot',
                            4: 'Label'},
                   inplace=True)
    print('in total: ' + str(len(df_data)) + ' drug-target pairs')

    df_data = encode_drug(df_data, drug_encoding)
    df_data = encode_protein(df_data, target_encoding)
    df_data = encode_structure(df_data, "0/1/2")
    df_data = encode_knot(df_data, "0/1/2")

    # dti split
    print('splitting dataset...')
    test = df_data.sample(frac=1, replace=False, random_state=random_seed)

    print('Done.')
    return test.reset_index(drop=True)


#########################################################
## 具体进行编码化的过程
#########################################################

# '?' padding
#amino_char = ['?', 'A', 'C', 'B', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'O','N', 'Q', 'P', 'S', 'R', 'U', 'T', 'W', 'V', 'Y', 'X', 'Z']  # 26
#nucleic_char = ['A', 'U', 'C', 'G', 'X', 'N']  # 7
#nucleic_char = ['A', 'U', 'C', 'G', 'N']  # 5
nucleic_char = ['A', 'U', 'C', 'G']  # 4

smiles_char = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '?', '#', '%', ')', '(', '[', ']', '+', '-', '.',
               '=', 'A', 'C', 'B', 'E', 'D', 'G', 'F', 'I',
               'H', 'K', 'M', 'L', 'O', 'N', 'P', 'S', 'R', 'U', 'T', 'W', 'V',
               'Y', '[', 'Z', ']', '_', 'a', 'c', 'b', 'e', 'd', 'g', 'f', 'i',
               'h', 'm', 'l', 'o', 'n', 's', 'r', 'u', 't', 'y']  # 65

#structure_char = ['.', '(', ')', '?']  # 5
structure_char = ['.', '(']  #2
#structure_char = ['.', '(', '[']  #3
#structure_char = ['.', '(', ')']  #3

knot_char = ['.', '[']  #3
#knot_char = [ '.', '[' ,'{','<','A','B','C','D','E']

from sklearn.preprocessing import OneHotEncoder

enc_protein = OneHotEncoder().fit(np.array(nucleic_char).reshape(-1, 1))
enc_drug = OneHotEncoder().fit(np.array(smiles_char).reshape(-1, 1))

MAX_SEQ_PROTEIN = 31
MAX_SEQ_DRUG = 100

def trans_protein(x):
    '''
    截取31bp长的RNA序列
    '''
    temp = list(x.upper())
    temp = [i if i in nucleic_char else 'N' for i in temp]
    if len(temp) < MAX_SEQ_PROTEIN:
        temp = temp + ['N'] * (MAX_SEQ_PROTEIN - len(temp))
    else:
        temp = temp[:MAX_SEQ_PROTEIN]
    return temp

def trans_structure(x):
    '''
    截取31bp长的RNA结构信息
    '''
    temp = list(x.upper())
    temp = [i if i in structure_char else '?' for i in temp]
    if len(temp) < MAX_SEQ_PROTEIN:
        temp = temp + ['?'] * (MAX_SEQ_PROTEIN - len(temp))
    else:
        temp = temp[:MAX_SEQ_PROTEIN]
    return temp

def trans_knot(x):
    '''
    截取31bp长的RNA knots
    '''
    temp = list(x.upper())
    temp = [i if i in knot_char else '?' for i in temp]
    if len(temp) < MAX_SEQ_PROTEIN:
        temp = temp + ['?'] * (MAX_SEQ_PROTEIN - len(temp))
    else:
        temp = temp[:MAX_SEQ_PROTEIN]
    return temp

def trans_drug(x):
    '''
    截取100原子单位的drug smile信息
    '''
    temp = list(x)
    temp = [i if i in smiles_char else '?' for i in temp]
    if len(temp) < MAX_SEQ_DRUG:
        temp = temp + ['?'] * (MAX_SEQ_DRUG - len(temp))
    else:
        temp = temp[:MAX_SEQ_DRUG]
    return temp


"""
def protein_2_embed(x):
    return enc_protein.transform(np.array(x).reshape(-1, 1)).toarray().T

def drug_2_embed(x):
    return enc_drug.transform(np.array(x).reshape(-1, 1)).toarray().T
"""


def save_dict(path, obj):
    with open(os.path.join(path, 'config.pkl'), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_dict(path):
    with open(os.path.join(path, 'config.pkl'), 'rb') as f:
        return pickle.load(f)


# ====================tranformer================================
import codecs
from subword_nmt.apply_bpe import BPE

# ESPF encoding
# vocab_path = f"smartnet/smartnet_classification/ESPF/drug_codes_chembl_freq_1500.txt"
vocab_path = f"smartnet/smartnet_classification/ESPF/drug_code.txt"
bpe_codes_drug = codecs.open(vocab_path)
# dbpe = BPE(bpe_codes_drug, merges=-1, separator='')
dbpe = BPE(bpe_codes_drug, merges=-1, separator='')

# sub_csv = pd.read_csv(f"smartnet/smartnet_classification/ESPF/subword_units_map_chembl_freq_1500.csv")
sub_csv = pd.read_csv(f"smartnet/smartnet_classification/ESPF/subword_units_drug.csv")
idx2word_d = sub_csv['index'].values
words2idx_d = dict(zip(idx2word_d, range(0, len(idx2word_d))))

# vocab_path = f"smartnet/smartnet_classification/ESPF/protein_codes_uniprot_2000.txt"
vocab_path = f"smartnet/smartnet_classification/ESPF/rna_code.txt"
bpe_codes_protein = codecs.open(vocab_path)
pbpe = BPE(bpe_codes_protein, merges=0, separator='')

# sub_csv = pd.read_csv(f"smartnet/smartnet_classification/ESPF/subword_units_map_uniprot_2000.csv")
sub_csv = pd.read_csv(f"smartnet/smartnet_classification/ESPF/subword_units_rna.csv")
idx2word_p = sub_csv['index'].values
words2idx_p = dict(zip(idx2word_p, range(0, len(idx2word_p))))


def drug2emb_encoder(x):
    # max_d = 50
    max_d = 100
    t1 = dbpe.process_line(x).split()  # split
    try:
        i1 = np.asarray([words2idx_d[i] for i in t1])  # index
    except:
        i1 = np.array([0])

    l = len(i1)

    if l < max_d:
        i = np.pad(i1, (0, max_d - l), 'constant', constant_values=0)
        input_mask = ([1] * l) + ([0] * (max_d - l))

    else:
        i = i1[:max_d]
        input_mask = [1] * max_d

    return i, np.asarray(input_mask)

def protein2emb_encoder(x):
    # max_p = 545
    max_p = 31
    t1 = pbpe.process_line(x).split()  # split
    try:
        i1 = np.asarray([words2idx_p[i] for i in t1])  # index
    except:
        i1 = np.array([0])

    l = len(i1)

    if l < max_p:
        i = np.pad(i1, (0, max_p - l), 'constant', constant_values=0)
        input_mask = ([1] * l) + ([0] * (max_p - l))
    else:
        i = i1[:max_p]
        input_mask = [1] * max_p

    return i, np.asarray(input_mask)

# =====================predicitoon==========================
def data_process3(X_drug=None, X_target=None, drug_encoding=None, target_encoding=None):
    y = [0] * int(len(X_drug) * len(X_target))  # create temp y for compatitibility

    if isinstance(X_target, str):
        X_target = [X_target]

    d_n = []
    t_n = []

    if len(X_target) == 1 and len(X_drug) > 1:
        # one target high throughput screening setting
        t_n = np.tile(X_target, (length_func(X_drug),))
        d_n = X_drug
    elif len(X_drug) == 1 and len(X_target) > 1:
        # one drug high throughput screening setting
        d_n = np.tile(X_drug, (length_func(X_target),))
        t_n = X_target

    elif len(X_drug) == 1 and len(X_target) == 1:
        d_n = X_drug
        t_n = X_target
    elif len(X_drug) > 1 and len(X_target) > 1:
        for i in range(len(X_drug)):
            for j in range(len(X_target)):
                d_n.append(X_drug[i])
                t_n.append(X_target[j])

    df_data = pd.DataFrame(zip(d_n, t_n, y))

    df_data.rename(columns={0: 'SMILES',
                            1: 'Target Sequence',
                            2: 'Label'},
                   inplace=True)
    print('in total: ' + str(len(df_data)) + ' drug-target pairs')

    df_data = encode_drug(df_data, drug_encoding)
    df_data = encode_protein(df_data, target_encoding)

    print('Done.')
    return df_data.reset_index(drop=True)

def data_process_repurpose_virtual_screening(drug, target, drug_encoding, target_encoding):
    df = data_process3(drug, target, drug_encoding=drug_encoding, target_encoding=target_encoding)

    return df

# =====================predicitoon_structure=================
def data_process3_ss(X_drug=None, X_target=None, X_structure=None,  drug_encoding=None, target_encoding=None):
    y = [0] * int(len(X_drug) * len(X_target))  # create temp y for compatitibility

    if isinstance(X_target, str):
        X_target = [X_target]

    d_n = []
    t_n = []
    s_n = []

    if len(X_target) == 1 and len(X_drug) > 1:
        # one target high throughput screening setting
        t_n = np.tile(X_target, (length_func(X_drug),))
        d_n = X_drug
        s_n = np.tile(X_structure, (length_func(X_drug),))
    elif len(X_drug) == 1 and len(X_target) > 1:
        # one drug high throughput screening setting
        d_n = np.tile(X_drug, (length_func(X_target),))
        t_n = X_target
        s_n = X_structure
    elif len(X_drug) == 1 and len(X_target) == 1:
        d_n = X_drug
        t_n = X_target
        s_n = X_structure
    elif len(X_drug) > 1 and len(X_target) > 1:
        for i in range(len(X_drug)):
            for j in range(len(X_target)):
                d_n.append(X_drug[i])
                t_n.append(X_target[j])
                s_n.append(X_structure[j])

    df_data = pd.DataFrame(zip(d_n, t_n, s_n, y))




    df_data.rename(columns={0: 'SMILES',
                            1: 'Target Sequence',
                            2: 'Target Structure',
                            3: 'Label'},
                   inplace=True)
    print('in total: ' + str(len(df_data)) + ' drug-target pairs')

    df_data = encode_drug(df_data, drug_encoding)
    df_data = encode_protein(df_data, target_encoding)
    df_data = encode_structure(df_data, "0/1/2")

    print('Done.')
    return df_data.reset_index(drop=True)

def data_process_repurpose_virtual_screening_ss(drug, target, structure, drug_encoding, target_encoding):
    df = data_process3_ss(drug, target, structure, drug_encoding=drug_encoding, target_encoding=target_encoding)
    return df

