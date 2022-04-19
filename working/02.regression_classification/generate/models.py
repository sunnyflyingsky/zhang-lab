import torch
from torch.autograd import Variable
import torch.nn.functional as F
from torch.utils import data
from torch.utils.data import SequentialSampler
from torch import nn

# from tqdm import tqdm
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from time import time
from sklearn.metrics import mean_squared_error, roc_auc_score, average_precision_score, f1_score, log_loss
# from lifelines.utils import concordance_index
from scipy.stats import pearsonr
import pickle

"""
import random
import os
seed=1
torch.set_num_threads(1)  # Suggested for issues with deadlocks, etc.
random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.enabled = True
"""

import copy
from prettytable import PrettyTable

import os

os.environ['CUDA_VISIBLE_DEVICES'] = '0'
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

cpu_num = 1  # 这里设置成你想运行的CPU个数
os.environ['OMP_NUM_THREADS'] = str(cpu_num)
os.environ['OPENBLAS_NUM_THREADS'] = str(cpu_num)
os.environ['MKL_NUM_THREADS'] = str(cpu_num)
os.environ['VECLIB_MAXIMUM_THREADS'] = str(cpu_num)
os.environ['NUMEXPR_NUM_THREADS'] = str(cpu_num)
torch.set_num_threads(cpu_num)

"""
class Classifier(nn.Sequential):
	def __init__(self, model_drug, model_protein, **config):
		super(Classifier, self).__init__()
		self.input_dim_drug = config['hidden_dim_drug']
		self.input_dim_protein = config['hidden_dim_protein']

		self.model_drug = model_drug
		self.model_protein = model_protein

		self.dropout = nn.Dropout(0.1)

		self.hidden_dims = config['cls_hidden_dims']
		layer_size = len(self.hidden_dims) + 1
		dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]

		self.predictor = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

	def forward(self, v_D, v_P):
		# each encoding
		v_D = self.model_drug(v_D)  # 128
		v_P = self.model_protein(v_P)  # 128
		# concatenate and classify
		v_f = torch.cat((v_D, v_P), 1)  # 256
		for i, l in enumerate(self.predictor):
			if i == (len(self.predictor) - 1):
				v_f = l(v_f)
			else:
				v_f = F.relu(self.dropout(l(v_f)))
		return v_f
"""


# ============================================CNN+CNN=============================================================

class CNN_CNN(nn.Sequential):

    def __init__(self, **config):
        super(CNN_CNN, self).__init__()

        ###drug
        in_ch_d = [63] + config['cnn_drug_filters']
        kernels_d = config['cnn_drug_kernels']
        layer_size_d = len(config['cnn_drug_filters'])
        self.convd = nn.ModuleList([nn.Conv1d(in_channels=in_ch_d[i],
                                              out_channels=in_ch_d[i + 1],
                                              kernel_size=kernels_d[i]) for i in range(layer_size_d)])
        self.convd = self.convd.double()

        n_size_d = self._get_conv_output_d((63, 100))
        self.fcd = nn.Linear(n_size_d, config['hidden_dim_drug'])
        # self.fcd = nn.Linear(in_ch_d[-1], config['hidden_dim_drug'])

        ###protein
        in_ch_t = [6] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])
        self.convt = self.convt.double()
        n_size_t = self._get_conv_output_t((6, 31))
        self.fct = nn.Linear(n_size_t, config['hidden_dim_protein'])
        # self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

    def _get_conv_output_d(self, shape):
        bs = 1
        input_d = Variable(torch.rand(bs, *shape))
        output_feat_d = self._forward_features_d(input_d.double())
        n_size_d = output_feat_d.data.view(bs, -1).size(1)
        return n_size_d

    def _forward_features_d(self, x):
        for l in self.convd:
            x = F.relu(l(x))
        x = F.adaptive_max_pool1d(x, output_size=1)
        return x

    def _get_conv_output_t(self, shape):
        bs = 1
        input_t = Variable(torch.rand(bs, *shape))
        output_feat_t = self._forward_features_t(input_t.double())
        n_size_t = output_feat_t.data.view(bs, -1).size(1)
        return n_size_t

    def _forward_features_t(self, x):
        for l in self.convt:
            x = F.relu(l(x))
        x = F.adaptive_max_pool1d(x, output_size=1)
        return x

    def forward(self, v_D, v_P):

        for l in self.convd:
            v_D = F.relu(l(v_D.double()))
        v_D = F.adaptive_max_pool1d(v_D, output_size=1)
        v_D = v_D.view(v_D.size(0), -1)
        v_D = self.fcd(v_D.float())

        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())

        # concatenate and classify

        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


# ====================================================CNN_RNN+RNN====================================================


class RNN_RNN(nn.Sequential):
    def __init__(self, **config):
        super(RNN_RNN, self).__init__()

        ###drug
        in_ch_d = [63] + config['cnn_drug_filters']
        self.in_ch_d = in_ch_d[-1]
        kernels_d = config['cnn_drug_kernels']
        layer_size_d = len(config['cnn_drug_filters'])
        self.convd = nn.ModuleList([nn.Conv1d(in_channels=in_ch_d[i],
                                              out_channels=in_ch_d[i + 1],
                                              kernel_size=kernels_d[i]) for i in range(layer_size_d)])
        self.convd = self.convd.double()
        n_size_d = self._get_conv_output_d((63, 100))  # auto get the seq_len of CNN output

        if config['rnn_Use_GRU_LSTM_drug'] == 'LSTM':
            self.rnnd = nn.LSTM(input_size=in_ch_d[-1],
                                hidden_size=config['rnn_drug_hid_dim'],
                                num_layers=config['rnn_drug_n_layers'],
                                batch_first=True,
                                bidirectional=config['rnn_drug_bidirectional'])

        elif config['rnn_Use_GRU_LSTM_drug'] == 'GRU':
            self.rnnd = nn.GRU(input_size=in_ch_d[-1],
                               hidden_size=config['rnn_drug_hid_dim'],
                               num_layers=config['rnn_drug_n_layers'],
                               batch_first=True,
                               bidirectional=config['rnn_drug_bidirectional'])
        else:
            raise AttributeError('Please use LSTM or GRU.')
        direction_d = 2 if config['rnn_drug_bidirectional'] else 1
        self.rnnd = self.rnnd.double()
        self.fcd = nn.Linear(config['rnn_drug_hid_dim'] * direction_d * n_size_d, config['hidden_dim_drug'])

        ###target
        in_ch_t = [6] + config['cnn_target_filters']
        self.in_ch_t = in_ch_t[-1]
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])
        self.convt = self.convt.double()
        n_size_t = self._get_conv_output_t((6, 31))

        if config['rnn_Use_GRU_LSTM_target'] == 'LSTM':
            self.rnnt = nn.LSTM(input_size=in_ch_t[-1],
                                hidden_size=config['rnn_target_hid_dim'],
                                num_layers=config['rnn_target_n_layers'],
                                batch_first=True,
                                bidirectional=config['rnn_target_bidirectional'])

        elif config['rnn_Use_GRU_LSTM_target'] == 'GRU':
            self.rnnt = nn.GRU(input_size=in_ch_t[-1],
                               hidden_size=config['rnn_target_hid_dim'],
                               num_layers=config['rnn_target_n_layers'],
                               batch_first=True,
                               bidirectional=config['rnn_target_bidirectional'])
        else:
            raise AttributeError('Please use LSTM or GRU.')
        direction_t = 2 if config['rnn_target_bidirectional'] else 1
        self.rnnt = self.rnnt.double()
        self.fct = nn.Linear(config['rnn_target_hid_dim'] * direction_t * n_size_t, config['hidden_dim_protein'])

        self.config = config

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

    def _get_conv_output_d(self, shape):
        bs = 1
        input_d = Variable(torch.rand(bs, *shape))
        output_feat_d = self._forward_features_d(input_d.double())
        n_size_d = output_feat_d.data.view(bs, self.in_ch_d, -1).size(2)
        return n_size_d

    def _forward_features_d(self, x):
        for l in self.convd:
            x = F.relu(l(x))
        return x

    def _get_conv_output_t(self, shape):
        bs = 1
        input_t = Variable(torch.rand(bs, *shape))
        output_feat_t = self._forward_features_t(input_t.double())
        n_size_t = output_feat_t.data.view(bs, self.in_ch_t, -1).size(2)
        return n_size_t

    def _forward_features_t(self, x):
        for l in self.convt:
            x = F.relu(l(x))
        return x

    def forward(self, v_D, v_P):
        for l in self.convd:
            v_D = F.relu(l(v_D.double()))
        batch_size_d = v_D.size(0)
        v_D = v_D.view(v_D.size(0), v_D.size(2), -1)
        if self.config['rnn_Use_GRU_LSTM_drug'] == 'LSTM':
            direction_d = 2 if self.config['rnn_drug_bidirectional'] else 1
            h0_d = torch.randn(self.config['rnn_drug_n_layers'] * direction_d, batch_size_d,
                               self.config['rnn_drug_hid_dim']).to(device)
            c0_d = torch.randn(self.config['rnn_drug_n_layers'] * direction_d, batch_size_d,
                               self.config['rnn_drug_hid_dim']).to(device)
            v_D, (hn_d, cn_d) = self.rnnd(v_D.double(), (h0_d.double(), c0_d.double()))
        else:
            # GRU
            direction_d = 2 if self.config['rnn_drug_bidirectional'] else 1
            h0_d = torch.randn(self.config['rnn_drug_n_layers'] * direction_d, batch_size_d,
                               self.config['rnn_drug_hid_dim']).to(device)
            v_D, hn_d = self.rnnd(v_D.double(), h0_d.double())

        v_D = torch.flatten(v_D, 1)
        v_D = self.fcd(v_D.float())

        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        batch_size_t = v_P.size(0)
        v_P = v_P.view(v_P.size(0), v_P.size(2), -1)
        if self.config['rnn_Use_GRU_LSTM_target'] == 'LSTM':
            direction_t = 2 if self.config['rnn_target_bidirectional'] else 1
            h0_t = torch.randn(self.config['rnn_target_n_layers'] * direction_t, batch_size_t,
                               self.config['rnn_target_hid_dim']).to(device)
            c0_t = torch.randn(self.config['rnn_target_n_layers'] * direction_t, batch_size_t,
                               self.config['rnn_target_hid_dim']).to(device)
            v_P, (hn_t, cn_t) = self.rnnt(v_P.double(), (h0_t.double(), c0_t.double()))
        else:
            # GRU
            direction_t = 2 if self.config['rnn_target_bidirectional'] else 1
            h0_t = torch.randn(self.config['rnn_target_n_layers'] * direction_t, batch_size_t,
                               self.config['rnn_target_hid_dim']).to(device)
            v_P, hn_t = self.rnnt(v_P.double(), h0_t.double())

        v_P = torch.flatten(v_P, 1)
        v_P = self.fct(v_P.float())

        # concatenate and classify

        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


# ====================================================CNN+GCN====================================================


class CNN_GCN(nn.Sequential):

    # def __init__(self, in_feats, hidden_feats=None, activation=None, predictor_dim=None, **config):
    def __init__(self, in_feats, **config):
        super(CNN_GCN, self).__init__()
        from dgllife.model.gnn.gcn import GCN
        from dgllife.model.readout.weighted_sum_and_max import WeightedSumAndMax

        # drug
        # in_feats=74
        predictor_dim = config['hidden_dim_drug']
        hidden_feats = [config['gnn_hid_dim_drug']] * config['gnn_num_layers']
        activation = [config['gnn_activation']] * config['gnn_num_layers']

        self.gnn = GCN(in_feats=in_feats,
                       hidden_feats=hidden_feats,
                       activation=activation
                       )
        gnn_out_feats = self.gnn.hidden_feats[-1]
        self.readout = WeightedSumAndMax(gnn_out_feats)
        self.transform = nn.Linear(self.gnn.hidden_feats[-1] * 2, predictor_dim)

        ###protein
        in_ch_t = [6] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])

        self.convt = self.convt.double()

        self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

    def forward(self, v_D, v_P):
        #import pdb;pdb.set_trace()

        v_D = v_D.to(device)

        #feats = v_D.ndata.pop('h')
        #node_feats = self.gnn(v_D, feats)

        feats= v_D.ndata['h']
        node_feats = self.gnn(v_D, feats)


        graph_feats = self.readout(v_D, node_feats)
        v_D = self.transform(graph_feats)

        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())

        # concatenate and classify
        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


# ====================================================CNN+GAT====================================================


class CNN_GAT(nn.Sequential):

    # def __init__(self, in_feats, hidden_feats=None, activation=None, predictor_dim=None, **config):
    def __init__(self, in_feats, **config):
        super(CNN_GAT, self).__init__()
        from dgllife.model.gnn.gat import GAT
        from dgllife.model.readout.weighted_sum_and_max import WeightedSumAndMax

        # drug
        # in_feats=74
        predictor_dim = config['hidden_dim_drug']
        hidden_feats = [config['gnn_hid_dim_drug']] * config['gnn_num_layers']
        activations = [config['gnn_activation']] * config['gnn_num_layers']
        num_heads = [config['gat_num_heads']] * config['gnn_num_layers']
        feat_drops = [config['gat_feat_drops']] * config['gnn_num_layers']
        attn_drops = [config['gat_attn_drops']] * config['gnn_num_layers']
        residuals = [True] * config['gnn_num_layers']
        agg_modes= ['flatten'] * config['gnn_num_layers']

        self.gnn = GAT(in_feats=in_feats,
                       hidden_feats=hidden_feats,
                       activations=activations,
                       num_heads=num_heads,
                       feat_drops=feat_drops,
                       attn_drops=attn_drops,
                       residuals=residuals,
                       agg_modes=agg_modes
                       )
        gnn_out_feats = self.gnn.hidden_feats[-1]* num_heads[-1]
        self.readout = WeightedSumAndMax(gnn_out_feats)

        self.transform = nn.Linear(self.gnn.hidden_feats[-1] * 2 * num_heads[-1], predictor_dim)

        ###protein
        in_ch_t = [6] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])

        self.convt = self.convt.double()

        self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

    def forward(self, v_D, v_P):

        v_D = v_D.to(device)
        feats = v_D.ndata.pop('h')
        node_feats = self.gnn(v_D, feats)
        # print(v_D)
        graph_feats = self.readout(v_D, node_feats)
        v_D = self.transform(graph_feats)

        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())

        # concatenate and classify
        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


# ====================================================CNN+GIN_attriMasking====================================================

class CNN_GIN_AttrMasking(nn.Module):
    ## adapted from https://github.com/awslabs/dgl-lifesci/blob/2fbf5fd6aca92675b709b6f1c3bc3c6ad5434e96/examples/property_prediction/moleculenet/utils.py#L76
    def __init__(self, **config):
        super(CNN_GIN_AttrMasking, self).__init__()
        from dgllife.model import load_pretrained
        from dgl.nn.pytorch.glob import AvgPooling

        predictor_dim = config['hidden_dim_drug']
        ## this is fixed hyperparameters as it is a pretrained model
        self.gnn = load_pretrained('gin_supervised_masking')

        self.readout = AvgPooling()
        self.transform = nn.Linear(300, predictor_dim)


        ###protein
        #in_ch_t = [6] + config['cnn_target_filters']
        in_ch_t = [7] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])

        self.convt = self.convt.double()

        self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])



    def forward(self, v_D, v_P):
        v_D = v_D.to(device)
        node_feats = [
            v_D.ndata.pop('atomic_number'),
            v_D.ndata.pop('chirality_type')
        ]
        edge_feats = [
            v_D.edata.pop('bond_type'),
            v_D.edata.pop('bond_direction_type')
        ]

        node_feats = self.gnn(v_D, node_feats, edge_feats)
        graph_feats = self.readout(v_D, node_feats)
        v_D = self.transform(graph_feats)



        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())

        # concatenate and classify
        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f



# ====================================================CNN+GIN_ContextPred====================================================


class CNN_GIN_ContextPred(nn.Module):
    ## adapted from https://github.com/awslabs/dgl-lifesci/blob/2fbf5fd6aca92675b709b6f1c3bc3c6ad5434e96/examples/property_prediction/moleculenet/utils.py#L76
    def __init__(self, **config):
        super(CNN_GIN_ContextPred, self).__init__()
        from dgllife.model import load_pretrained
        from dgl.nn.pytorch.glob import AvgPooling


        predictor_dim = config['hidden_dim_drug']
        ## this is fixed hyperparameters as it is a pretrained model
        self.gnn = load_pretrained('gin_supervised_contextpred')

        self.readout = AvgPooling()
        self.transform = nn.Linear(300, predictor_dim)




        ###protein
        in_ch_t = [6] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])

        self.convt = self.convt.double()

        self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])


    def forward(self, v_D, v_P):
        v_D = v_D.to(device)
        node_feats = [
            v_D.ndata.pop('atomic_number'),
            v_D.ndata.pop('chirality_type')
        ]
        edge_feats = [
            v_D.edata.pop('bond_type'),
            v_D.edata.pop('bond_direction_type')
        ]

        node_feats = self.gnn(v_D, node_feats, edge_feats)
        graph_feats = self.readout(v_D, node_feats)
        v_D = self.transform(graph_feats)


        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())

        # concatenate and classify
        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f
# ====================================================CNN+NeuralFP====================================================

class CNN_NeuralFP(nn.Sequential):
    # def __init__(self, in_feats, hidden_feats=None, activation=None, predictor_dim=None, **config):
    def __init__(self, in_feats, **config):
        super(CNN_NeuralFP, self).__init__()
        from dgllife.model.gnn.nf import NFGNN
        from dgllife.model.readout.sum_and_max import SumAndMax

        # drug
        # in_feats=74

        hidden_feats = [config['gnn_hid_dim_drug']] * config['gnn_num_layers']
        max_degree = config['neuralfp_max_degree']
        activation = [config['gnn_activation']] * config['gnn_num_layers']
        predictor_hidden_size = config['neuralfp_predictor_hid_dim']
        predictor_dim = config['hidden_dim_drug']
        predictor_activation = config['neuralfp_predictor_activation']

        self.gnn = NFGNN(in_feats=in_feats,
                         hidden_feats=hidden_feats,
                         max_degree=max_degree,
                         activation=activation
                         )
        gnn_out_feats = self.gnn.gnn_layers[-1].out_feats
        self.node_to_graph = nn.Linear(gnn_out_feats, predictor_hidden_size)
        self.predictor_activation = predictor_activation

        self.readout = SumAndMax()
        self.transform = nn.Linear(predictor_hidden_size * 2, predictor_dim)

        ###protein
        in_ch_t = [6] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])

        self.convt = self.convt.double()

        self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

    def forward(self, v_D, v_P):

        v_D = v_D.to(device)
        # v_D = v_D
        feats = v_D.ndata.pop('h')
        node_feats = self.gnn(v_D, feats)
        node_feats = self.node_to_graph(node_feats)
        graph_feats = self.readout(v_D, node_feats)
        graph_feats = self.predictor_activation(graph_feats)
        v_D = self.transform(graph_feats)

        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())

        # concatenate and classify
        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


# ====================================================CNN+AttentiveFP====================================================


class CNN_AttentiveFP(nn.Sequential):

    # def __init__(self, in_feats, hidden_feats=None, activation=None, predictor_dim=None, **config):
    def __init__(self, **config):
        super(CNN_AttentiveFP, self).__init__()
        from dgllife.model.gnn import AttentiveFPGNN
        from dgllife.model.readout import AttentiveFPReadout
        # drug
        # in_feats=74

        node_feat_size = 39
        edge_feat_size = 11
        num_layers = config['gnn_num_layers']
        num_timesteps = config['attentivefp_num_timesteps']
        graph_feat_size = config['gnn_hid_dim_drug']
        predictor_dim = config['hidden_dim_drug']

        self.gnn = AttentiveFPGNN(node_feat_size=node_feat_size,
                                  edge_feat_size=edge_feat_size,
                                  num_layers=num_layers,
                                  graph_feat_size=graph_feat_size)

        self.readout = AttentiveFPReadout(feat_size=graph_feat_size,
                                          num_timesteps=num_timesteps)

        self.transform = nn.Linear(graph_feat_size, predictor_dim)

        ###protein
        in_ch_t = [5] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']
        layer_size_t = len(config['cnn_target_filters'])
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])

        self.convt = self.convt.double()

        self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

    def forward(self, v_D, v_P):

        v_D = v_D.to(device)
        # v_D = v_D
        node_feats = v_D.ndata.pop('h')
        edge_feats = v_D.edata.pop('e')
        node_feats = self.gnn(v_D, node_feats, edge_feats)
        graph_feats = self.readout(v_D, node_feats, False)
        v_D = self.transform(graph_feats)

        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())

        # concatenate and classify
        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


# ====================================================Transformer====================================================


from model_helper import *


class Trans_Trans(nn.Sequential):

    def __init__(self, **config):
        super(Trans_Trans, self).__init__()

        # self.emb_d = Embeddings(config['input_dim_drug'], config['transformer_emb_size_drug'], 50,config['transformer_dropout_rate'])

        self.emb_d = Embeddings(config['input_dim_drug'], config['transformer_emb_size_drug'], 100,
                                config['transformer_dropout_rate'])

        self.encoder_d = Encoder_MultipleLayers(config['transformer_n_layer_drug'],
                                                config['transformer_emb_size_drug'],
                                                config['transformer_intermediate_size_drug'],
                                                config['transformer_num_attention_heads_drug'],
                                                config['transformer_attention_probs_dropout'],
                                                config['transformer_hidden_dropout_rate'])

        self.emb_p = Embeddings(config['input_dim_protein'], config['transformer_emb_size_target'], 31,
                                config['transformer_dropout_rate'])
        self.encoder_p = Encoder_MultipleLayers(config['transformer_n_layer_target'],
                                                config['transformer_emb_size_target'],
                                                config['transformer_intermediate_size_target'],
                                                config['transformer_num_attention_heads_target'],
                                                config['transformer_attention_probs_dropout'],
                                                config['transformer_hidden_dropout_rate'])

        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])

    ### parameter v (tuple of length 2) is from utils.drug2emb_encoder
    def forward(self, v_D, v_P):
        e_d = v_D[0].long().to(device)
        e_mask_d = v_D[1].long().to(device)
        ex_e_mask_d = e_mask_d.unsqueeze(1).unsqueeze(2)
        ex_e_mask_d = (1.0 - ex_e_mask_d) * -10000.0
        emb_d = self.emb_d(e_d)
        encoded_layers_d = self.encoder_d(emb_d.float(), ex_e_mask_d.float())
        v_D = encoded_layers_d[:, 0]

        e_t = v_P[0].long().to(device)
        e_mask_t = v_P[1].long().to(device)
        ex_e_mask_t = e_mask_t.unsqueeze(1).unsqueeze(2)
        ex_e_mask_t = (1.0 - ex_e_mask_t) * -10000.0
        emb_t = self.emb_p(e_t)
        encoded_layers_t = self.encoder_p(emb_t.float(), ex_e_mask_t.float())
        v_P = encoded_layers_t[:, 0]

        # concatenate and classify

        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


# ====================================================CNN_CNN_inter===================================================


class Conv1dReLU(nn.Module):
    '''
    kernel_size=3, stride=1, padding=1
    kernel_size=5, stride=1, padding=2
    kernel_size=7, stride=1, padding=3
    '''

    def __init__(self, in_channels, out_channels, kernel_size, stride=1, padding=0):
        super().__init__()
        self.inc = nn.Sequential(
            nn.Conv1d(in_channels=in_channels, out_channels=out_channels, kernel_size=kernel_size, stride=stride,
                      padding=padding),
            nn.ReLU()
        )

    def forward(self, x):
        return self.inc(x)


class LinearReLU(nn.Module):
    def __init__(self, in_features, out_features, bias=True):
        super().__init__()
        self.inc = nn.Sequential(
            nn.Linear(in_features=in_features, out_features=out_features, bias=bias),
            nn.ReLU()
        )

    def forward(self, x):
        return self.inc(x)


class BilinearPooling(nn.Module):
    def __init__(self, in_channels, out_channels, c_m, c_n):
        super().__init__()

        self.convA = nn.Conv1d(in_channels, c_m, kernel_size=1, stride=1, padding=0)
        self.convB = nn.Conv1d(in_channels, c_n, kernel_size=1, stride=1, padding=0)
        self.linear = nn.Linear(c_m, out_channels, bias=True)


class MutualAttentation(nn.Module):
    def __init__(self, in_channels, att_size, c_m, c_n):
        super().__init__()
        self.bipool = BilinearPooling(in_channels, in_channels, c_m, c_n)
        self.linearS = nn.Linear(in_channels, att_size)
        self.linearT = nn.Linear(in_channels, att_size)


class CNN_CNN_inter(nn.Sequential):

    def __init__(self, **config):
        super(CNN_CNN_inter, self).__init__()

        in_ch_d = [63] + config['cnn_drug_filters']
        kernels_d = config['cnn_drug_kernels']

        in_ch_t = [6] + config['cnn_target_filters']
        kernels_t = config['cnn_target_kernels']

        self.prot_conv1 = Conv1dReLU(in_ch_t[0], in_ch_t[1], kernels_t[0])
        self.prot_conv2 = Conv1dReLU(in_ch_t[1], in_ch_t[2], kernels_t[1])
        self.prot_conv3 = Conv1dReLU(in_ch_t[2], in_ch_t[3], kernels_t[2])
        self.prot_pool = nn.AdaptiveMaxPool1d(1)

        self.drug_conv1 = Conv1dReLU(in_ch_d[0], in_ch_d[1], kernels_d[0])
        self.drug_conv2 = Conv1dReLU(in_ch_d[1], in_ch_d[2], kernels_d[1])
        self.drug_conv3 = Conv1dReLU(in_ch_d[2], in_ch_d[3], kernels_d[2])
        self.drug_pool = nn.AdaptiveMaxPool1d(1)

        self.prot_mut_att1 = MutualAttentation(in_ch_t[1], in_ch_t[1], in_ch_t[1], 8)
        self.drug_mut_att1 = MutualAttentation(in_ch_d[1], in_ch_d[1], in_ch_d[1], 8)

        self.prot_mut_att2 = MutualAttentation(in_ch_t[2], in_ch_t[1], in_ch_t[1], 8)
        self.drug_mut_att2 = MutualAttentation(in_ch_d[2], in_ch_d[1], in_ch_d[1], 8)

        self.prot_mut_att3 = MutualAttentation(in_ch_t[3], in_ch_t[1], in_ch_t[1], 8)
        self.drug_mut_att3 = MutualAttentation(in_ch_d[3], in_ch_d[1], in_ch_d[1], 8)

        self.linear1 = LinearReLU(in_ch_d[3] + in_ch_t[3], 1024)
        self.drop1 = nn.Dropout(0.1)
        self.linear2 = LinearReLU(1024, 1024)
        self.drop2 = nn.Dropout(0.1)
        self.linear3 = LinearReLU(1024, 512)
        self.drop3 = nn.Dropout(0.1)
        self.out_layer = nn.Linear(512, 1)

    def forward(self, v_D, v_P):
        prot_x = self.prot_conv1(v_P)
        drug_x = self.drug_conv1(v_D)
        prot_x_g = self.prot_mut_att1(drug_x, prot_x)
        drug_x_g = self.drug_mut_att1(prot_x, drug_x)

        prot_x = self.prot_conv2(prot_x_g)
        drug_x = self.drug_conv2(drug_x_g)
        prot_x_g = self.prot_mut_att2(drug_x, prot_x)
        drug_x_g = self.drug_mut_att2(prot_x, drug_x)

        prot_x = self.prot_conv3(prot_x_g)
        drug_x = self.drug_conv3(drug_x_g)
        prot_x_g = self.prot_mut_att3(drug_x, prot_x)
        drug_x_g = self.drug_mut_att3(prot_x, drug_x)

        prot_x = self.prot_pool(prot_x_g).squeeze(-1)
        drug_x = self.drug_pool(drug_x_g).squeeze(-1)

        x = torch.cat([prot_x, drug_x], dim=-1)
        x = self.linear1(x)
        x = self.drop1(x)
        x = self.linear2(x)
        x = self.drop2(x)
        x = self.linear3(x)
        x = self.drop3(x)
        x = self.out_layer(x)

        return x






# ====================================================Prism+GIN===================================================
class Conv2d(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, stride=1, relu=True, same_padding=False, bn=False):
        super(Conv2d, self).__init__()
        p0 = int((kernel_size[0] - 1) / 2) if same_padding else 0
        p1 = int((kernel_size[1] - 1) / 2) if same_padding else 0
        padding = (p0, p1)
        self.conv = nn.Conv2d(in_channels, out_channels, kernel_size, stride, padding=padding)
        self.bn = nn.BatchNorm2d(out_channels) if bn else None
        self.relu = nn.ReLU(inplace=True) if relu else None

    def forward(self, x):
        x = self.conv(x)
        if self.bn is not None:
            x = self.bn(x)
        if self.relu is not None:
            x = self.relu(x)
        return x

class SEBlock(nn.Module):
    def __init__(self, channel, reduction=2):
        super(SEBlock, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool2d(1)
        self.fc = nn.Sequential(
                nn.Linear(channel, channel // reduction),
                nn.ReLU(inplace=True),
                nn.Linear(channel // reduction, channel),
                nn.Sigmoid()
        )

    def forward(self, x):
        b, c, _, _ = x.size()
        y = self.avg_pool(x).view(b, c)
        y = self.fc(y).view(b, c, 1, 1)
        return y


class ResidualBlock1D(nn.Module):

    def __init__(self, planes, downsample=True):
        super(ResidualBlock1D, self).__init__()
        self.c1 = nn.Conv1d(planes,   planes,   kernel_size=1, stride=1, bias=False)
        self.b1 = nn.BatchNorm1d(planes)
        self.c2 = nn.Conv1d(planes,   planes*2, kernel_size=11, stride=1,padding=5, bias=False)
        #self.c2 = nn.Conv1d(planes,   planes*2, kernel_size=7, stride=1,padding=3, bias=False)
        self.b2 = nn.BatchNorm1d(planes*2)
        self.c3 = nn.Conv1d(planes*2, planes*8, kernel_size=1, stride=1, bias=False)
        self.b3 = nn.BatchNorm1d(planes * 8)
        self.downsample = nn.Sequential(
            nn.Conv1d(planes,   planes*8,   kernel_size=1, stride=1, bias=False),
            nn.BatchNorm1d(planes*8),
        )
        self.relu  = nn.ReLU(inplace=True)

    def forward(self, x):
        identity = x

        out = self.c1(x)
        out = self.b1(out)
        out = self.relu(out)

        out = self.c2(out)
        out = self.b2(out)
        out = self.relu(out)

        out = self.c3(out)
        out = self.b3(out)

        if self.downsample:
            identity = self.downsample(x)

        out += identity
        out = self.relu(out)

        return out

class ResidualBlock2D(nn.Module):

    def __init__(self, planes, kernel_size=(7,4), padding=(3,1), downsample=True):
        super(ResidualBlock2D, self).__init__()
        self.c1 = nn.Conv2d(planes,   planes,   kernel_size=1, stride=1, bias=False)
        self.b1 = nn.BatchNorm2d(planes)
        self.c2 = nn.Conv2d(planes,   planes*2, kernel_size=kernel_size, stride=1,
                     padding=padding, bias=False)
        self.b2 = nn.BatchNorm2d(planes*2)
        self.c3 = nn.Conv2d(planes*2, planes*4, kernel_size=1, stride=1, bias=False)
        self.b3 = nn.BatchNorm2d(planes * 4)
        self.downsample = nn.Sequential(
            nn.Conv2d(planes,   planes*4,   kernel_size=1, stride=1, bias=False),
            nn.BatchNorm2d(planes*4),
        )
        self.relu  = nn.ReLU(inplace=True)

    def forward(self, x):
        identity = x

        out = self.c1(x)
        out = self.b1(out)
        out = self.relu(out)

        out = self.c2(out)
        out = self.b2(out)
        out = self.relu(out)

        out = self.c3(out)
        out = self.b3(out)

        if self.downsample:
            identity = self.downsample(x)
        #print(out.shape,identity.shape)
        out += identity
        out = self.relu(out)

        return out
        
class Prism_GIN_AttrMasking(nn.Module):

    def __init__(self, **config):
        super(Prism_GIN_AttrMasking, self).__init__()
        from dgllife.model import load_pretrained
        from dgl.nn.pytorch.glob import AvgPooling

        predictor_dim = config['hidden_dim_drug']
        ## this is fixed hyperparameters as it is a pretrained model
        self.gnn = load_pretrained('gin_supervised_masking')

        self.readout = AvgPooling()
        self.transform = nn.Linear(300, predictor_dim)


        ###protein

        #in_ch_t = [6] + config['cnn_target_filters']
        #kernels_t = config['cnn_target_kernels']
        #layer_size_t = len(config['cnn_target_filters'])


        #h_p, h_k = 0, 1
        #self.n_features = 1

        #h_p, h_k = 1, 2
        #self.n_features = 2

        #h_p, h_k = 1, 3
        #self.n_features = 3

        #h_p, h_k = 1, 3
        #self.n_features = 4


        #h_p, h_k = 2, 5
        #self.n_features = 5

        h_p, h_k = 2, 5
        self.n_features = 6


        #h_p, h_k = 4, 9
        #self.n_features = 10


        #5-2,3-1 for 15bp.  15-k+2p+1=
        #kernal=3
        #pad=1


        #9-4,7-3,5-2,3-1 for 31bp
        kernal=7
        pad=3


        #11-5, 9-4, 7-3, 5-2,3-1 for 51bp.  51-k+2p+1=
        #kernal=11
        #pad=5



        base_channel = 8
        """
        self.convt = nn.ModuleList([nn.Conv1d(in_channels=in_ch_t[i],
                                              out_channels=in_ch_t[i + 1],
                                              kernel_size=kernels_t[i]) for i in range(layer_size_t)])



        self.convt = self.convt.double()

        self.fct = nn.Linear(in_ch_t[-1], config['hidden_dim_protein'])
        """
        """
        import random
        import os
        seed=1
        torch.set_num_threads(1)  # Suggested for issues with deadlocks, etc.
        random.seed(seed)
        os.environ['PYTHONHASHSEED'] = str(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.enabled = True
        """



        self.conv    = Conv2d(1, base_channel, kernel_size=(kernal,h_k), bn = True, same_padding=True)
        self.se      = SEBlock(base_channel)
        self.res2d   = ResidualBlock2D(base_channel, kernel_size=(kernal, h_k), padding=(pad, h_p))  #
        self.res1d   = ResidualBlock1D(base_channel*4) 
        self.avgpool = nn.AvgPool2d((1,self.n_features))
        self.gpool   = nn.AdaptiveAvgPool1d(1)
        #self.fct      = nn.Linear(base_channel*4*8, 1)
        self.fct      = nn.Linear(base_channel*4*8, config['hidden_dim_protein'])



        ###dropout
        self.dropout = nn.Dropout(0.1)

        ###Fully-connected
        self.input_dim_drug = config['hidden_dim_drug']
        self.input_dim_protein = config['hidden_dim_protein']
        self.hidden_dims = config['cls_hidden_dims']
        dims = [self.input_dim_drug + self.input_dim_protein] + self.hidden_dims + [1]  # [256,1024,1024,512,1]
        layer_size = len(self.hidden_dims) + 1

        self.fc = nn.ModuleList([nn.Linear(dims[i], dims[i + 1]) for i in range(layer_size)])



    def forward(self, v_D, v_P):
        v_D = v_D.to(device)
        node_feats = [
            v_D.ndata.pop('atomic_number'),
            v_D.ndata.pop('chirality_type')
        ]
        edge_feats = [
            v_D.edata.pop('bond_type'),
            v_D.edata.pop('bond_direction_type')
        ]

        node_feats = self.gnn(v_D, node_feats, edge_feats)
        graph_feats = self.readout(v_D, node_feats)
        v_D = self.transform(graph_feats)


        """
        for l in self.convt:
            v_P = F.relu(l(v_P.double()))
        v_P = F.adaptive_max_pool1d(v_P, output_size=1)
        v_P = v_P.view(v_P.size(0), -1)
        v_P = self.fct(v_P.float())
        """

        #print("input: "+str(v_P.shape))
        #print(v_P) #AUG
        x = self.conv(v_P)
        #print("conv(x): "+str(x.shape))
        #print(x) 

        x = F.dropout(x, 0.1, training=self.training)


        #print("conv_dropout(x): "+str(x.shape))
        z = self.se(x)
        #print("se(z): "+str(z.shape))
        #z= z.repeat((1,1,31,5))
        #print("se(z): "+str(z.shape))
        x = self.res2d(z*x)
        x = F.dropout(x, 0.5, training=self.training)
        x = self.avgpool(x)
        x = x.view(x.shape[0], x.shape[1], x.shape[2])
        x = self.res1d(x)
        x = F.dropout(x, 0.3, training=self.training)
        x = self.gpool(x)
        x = x.view(x.shape[0], x.shape[1])
        v_P = self.fct(x)



        # concatenate and classify
        v_f = torch.cat((v_D, v_P), 1)  # 256

        for i, l in enumerate(self.fc):
            if i == (len(self.fc) - 1):
                v_f = l(v_f)
            else:
                v_f = F.relu(self.dropout(l(v_f)))
        return v_f


