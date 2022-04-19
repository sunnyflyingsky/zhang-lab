import torch
from torch.autograd import Variable
import torch.nn.functional as F
from torch.utils import data
from torch.utils.data import SequentialSampler
from torch import nn

#from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from time import time
from sklearn.metrics import mean_squared_error, roc_auc_score, average_precision_score, f1_score, log_loss, confusion_matrix
#from lifelines.utils import concordance_index
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

from loader import *
#from SmartNet.model_helper import Encoder_MultipleLayers, Embeddings
from models import *

from torch.utils.tensorboard import SummaryWriter


from processor import data_process_repurpose_virtual_screening
from processor import data_process_repurpose_virtual_screening_ss


def save_dict(path, obj):
    with open(os.path.join(path, 'config.pkl'), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_dict(path):
	with open(os.path.join(path, 'config.pkl'), 'rb') as f:
		return pickle.load(f)

def initialize(**config):
    model = SmartNet(**config)
    return model

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
         lw=lw, label= method_name + ' (area = %0.3f)' % roc_auc[0])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    fontsize = 14
    plt.xlabel('False Positive Rate', fontsize = fontsize)
    plt.ylabel('True Positive Rate', fontsize = fontsize)
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
    #plt.plot([0,1], [no_skill, no_skill], linestyle='--')
    plt.plot(lr_recall, lr_precision, lw = 2, label= method_name + ' (area = %0.3f)' % average_precision_score(y_label, y_pred))
    fontsize = 14
    plt.xlabel('Recall', fontsize = fontsize)
    plt.ylabel('Precision', fontsize = fontsize)
    plt.title('Precision Recall Curve')
    plt.legend()
    plt.savefig(figure_file)
    return

# ================================================================
def dgl_collate_func(x):
    d, p, y = zip(*x)
    import dgl
    d = dgl.batch(d)
    return d, torch.tensor(p), torch.tensor(y)
# ================================================================

class SmartNet:

    def __init__(self, **config):


        import os
        os.environ['CUDA_VISIBLE_DEVICES'] = '0'
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        cpu_num = 1  # 这里设置成你想运行的CPU个数
        os.environ['OMP_NUM_THREADS'] = str(cpu_num)
        os.environ['OPENBLAS_NUM_THREADS'] = str(cpu_num)
        os.environ['MKL_NUM_THREADS'] = str(cpu_num)
        os.environ['VECLIB_MAXIMUM_THREADS'] = str(cpu_num)
        os.environ['NUMEXPR_NUM_THREADS'] = str(cpu_num)
        torch.set_num_threads(cpu_num)


        drug_encoding = config['drug_encoding']
        target_encoding = config['target_encoding']
        inter = config['inter']


        #self.model = Classifier(self.model_drug, self.model_protein, **config)


        #self.model = Classifier_GCN(**config)

        if (drug_encoding == "CNN" and target_encoding == "CNN"):
            #self.model = CNN_CNN(**config)
            self.model = CNN_CNN_mul_DTI(**config)
            #self.model = CNN_CNN_mul(**config)
        elif (drug_encoding == "RNN" and target_encoding == "RNN"):
            self.model = RNN_RNN(**config)
        elif (drug_encoding == "DGL_GCN" and target_encoding == "CNN"):
            self.model = CNN_GCN(**config)
        elif (drug_encoding == "DGL_GAT" and target_encoding == "CNN"):
            self.model = CNN_GAT(**config)
        elif (drug_encoding == "DGL_NeuralFP" and target_encoding == "CNN"):
            self.model = CNN_NeuralFP(**config)
        elif (drug_encoding == "DGL_AttentiveFP" and target_encoding == "CNN"):
            self.model = CNN_AttentiveFP(**config)
        elif (drug_encoding == "DGL_GIN_AttrMasking" and target_encoding == "CNN"):
            #self.model = CNN_GIN_AttrMasking(**config)
            self.model = CNN_GIN_AttrMasking_mul(**config)
        elif (drug_encoding == "DGL_GIN_AttrMasking" and target_encoding == "Prism"):
            #self.model = Prism_GIN_AttrMasking(**config)
            #self.model = Prism_GIN_AttrMasking_mul(**config)
            self.model = Prism_GIN_AttrMasking_TF(**config)
        elif (drug_encoding == "DGL_GIN_ContextPred" and target_encoding == "CNN"):
            self.model = CNN_GIN_ContextPred(**config)
        elif (drug_encoding == "Transformer" and target_encoding == "Transformer"):
            self.model = Trans_Trans(**config)
        elif (drug_encoding == "CNN" and target_encoding == "CNN" and inter==1):
            self.model = CNN_CNN_inter(**config)

        self.config = config

        self.drug_encoding = drug_encoding
        self.target_encoding = target_encoding
        self.result_folder = config['result_folder']
        if not os.path.exists(self.result_folder):
            os.mkdir(self.result_folder)

        if 'num_workers' not in self.config.keys():
            self.config['num_workers'] = 0
        if 'decay' not in self.config.keys():
            self.config['decay'] = 0

    def train(self, train, val, test):

        lr = self.config['LR']
        decay = self.config['decay']
        BATCH_SIZE = self.config['batch_size']
        train_epoch = self.config['train_epoch']
        loss_history = []

        self.model = self.model.to(self.device)

        #============================= Support multiple GPUs===========================================
        if torch.cuda.device_count() > 1:
            print("Let's use " + str(torch.cuda.device_count()) + " GPUs!")
            self.model = nn.DataParallel(self.model, dim=0)
        elif torch.cuda.device_count() == 1:
            print("Let's use " + str(torch.cuda.device_count()) + " GPU!")
        else:
            print("Let's use CPU/s!")

        #============================= Support multiple optimizers with parameters======================
        def adjust_learning_rate(init_lr, optimizer, iteration):
            lr = max(init_lr * (0.9 ** (iteration//5)), 0.0001)
            for param_group in optimizer.param_groups:
                param_group["lr"] = lr
            return lr
        
        opt = torch.optim.Adam(self.model.parameters(), lr=lr, weight_decay=decay)

        #============================= Mini Batch training==================================================
        print('--- Data Preparation ---')

        params = {'batch_size': BATCH_SIZE,
                  'shuffle': True,
                  'num_workers': self.config['num_workers'],
                  'drop_last': False}

        if self.drug_encoding in ['DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking', 'DGL_GIN_ContextPred', 'DGL_AttentiveFP']:
            params['collate_fn'] = dgl_collate_func

        training_generator = data.DataLoader(data_process_loader(train.index.values, train.Label.values, train, **self.config), **params)
        validation_generator = data.DataLoader(data_process_loader(val.index.values, val.Label.values, val, **self.config), **params)
        params_test = {'batch_size': BATCH_SIZE,
                       'shuffle': False,
                       'num_workers': self.config['num_workers'],
                       'drop_last': False,
                       'sampler': SequentialSampler(data_process_loader(test.index.values, test.Label.values, test, **self.config))}

        if self.drug_encoding in ['DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking', 'DGL_GIN_ContextPred','DGL_AttentiveFP']:
            params_test['collate_fn'] = dgl_collate_func
        testing_generator = data.DataLoader(data_process_loader(test.index.values, test.Label.values, test, **self.config), **params_test)

        # early stopping
        max_auc = 0
        model_max = copy.deepcopy(self.model)

        valid_metric_record = []
        valid_metric_header = ["# epoch"]
        valid_metric_header.extend(["AUROC", "AUPRC", "F1"])
        table = PrettyTable(valid_metric_header)
        float2str = lambda x: '%0.4f' % x

        print('--- Go for Training ---')
        writer = SummaryWriter()
        t_start = time()
        iteration_loss = 0

        for epo in range(train_epoch):
            #lr=adjust_learning_rate(0.001, opt, epo)
            print(lr, epo)
            for i, (v_d, v_p, label) in enumerate(training_generator):
                #print(v_p.shape) #torch.Size([64, 5, 31])
                #print(type(v_p))
                #print(v_p)
                if self.target_encoding == 'Transformer':
                    v_p = v_p
                else:
                    v_p = v_p.float().to(self.device)
                    #v_p = v_p.float().to(self.device)
                if self.drug_encoding in ['Transformer', 'DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking','DGL_GIN_ContextPred', 'DGL_AttentiveFP']:
                    v_d = v_d
                else:
                    v_d = v_d.float().to(self.device)

                score = self.model(v_d, v_p)
                #print(label)
                #print(score)
                label = Variable(torch.from_numpy(np.array(label)).float()).to(self.device)

                # ===========loss function============

                loss_fct = torch.nn.BCELoss()
                m = torch.nn.Sigmoid()
                n = torch.squeeze(m(score), 1) ##去掉维数为1的维度
                loss = loss_fct(n, label)
                # ====================================


                # ===========loss function2============
                #loss_fct = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor(2)) #
                #m = torch.nn.Sigmoid()
                #n = torch.squeeze(score, 1)
                #loss = loss_fct(n, label)
                # ====================================
                loss_history.append(loss.item())
                writer.add_scalar("Loss/train", loss.item(), iteration_loss)
                iteration_loss += 1

                opt.zero_grad() #所有参数梯度为0
                loss.backward() #误差反向传递，计算梯度
                opt.step() #以学习效率优化梯度

                if (i % 100 == 0):
                    t_now = time()
                    print('Training at Epoch ' + str(epo + 1) + ' iteration ' + str(i) + \
                          ' with loss ' + str(loss.cpu().detach().numpy())[:7] + \
                          ". Total time " + str(int(t_now - t_start) / 3600)[:7] + " hours")
                    ### record total run time

            ##### validate, select the best model up to now
            with torch.set_grad_enabled(False):
                ## binary: ROC-AUC, PR-AUC, F1, cross-entropy loss
                y_pred = []
                y_label = []
                loss_lst=[]
                self.model.eval()
                for i, (v_d, v_p, label) in enumerate(validation_generator):
                    if self.drug_encoding in ['Transformer', "MPNN", 'Transformer', 'DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking','DGL_GIN_ContextPred', 'DGL_AttentiveFP']:
                        v_d = v_d
                    else:
                        v_d = v_d.float().to(self.device)
                    if self.target_encoding == 'Transformer':
                        v_p = v_p
                    else:
                        v_p = v_p.float().to(self.device)
                    score = self.model(v_d, v_p)
                    label = Variable(torch.from_numpy(np.array(label)).float()).to(self.device)
                    # binary
                    #m = torch.nn.Sigmoid()

                    # ===========loss function1============
                    loss_fct = torch.nn.BCELoss()
                    m = torch.nn.Sigmoid()
                    n = torch.squeeze(m(score), 1)
                    loss = loss_fct(n, label)
                    # ====================================

                    # ===========loss function2============
                    #loss_fct = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor(2)) #
                    #m = torch.nn.Sigmoid()
                    #n = torch.squeeze(score, 1)
                    #loss = loss_fct(n, label)
                    # ====================================

                    loss_lst.append(loss.item())

                    logits = torch.squeeze(m(score)).detach().cpu().numpy()
                    label_ids = label.to('cpu').numpy()
                    y_label = y_label + label_ids.flatten().tolist()
                    y_pred = y_pred + logits.flatten().tolist()
                    outputs = np.asarray([1 if i else 0 for i in (np.asarray(y_pred) >= 0.5)])
                self.model.train()

                ## Output
                auc = roc_auc_score(y_label, y_pred)
                auprc = average_precision_score(y_label, y_pred)
                f1 = f1_score(y_label, outputs)
                #loss = log_loss(y_label, outputs)
                loss_lst = np.array(loss_lst)

                lst = ["epoch " + str(epo)] + list(map(float2str, [auc, auprc, f1]))
                valid_metric_record.append(lst)
                if auc > max_auc:
                    model_max = copy.deepcopy(self.model)  # # load early stopped model
                    max_auc = auc
                    best_epo = epo
                if epo - best_epo > 20:
                    print("Early stop at %d"%(epo + 1))
                    break
                if  epo == best_epo:
                    print('Validation at Epoch ' + str(epo + 1) + ', AUROC: ' + str(auc)[:7] + ' , AUPRC: ' + str(auprc)[:7] + ' , F1: ' + str(f1)[:7] + ' , Cross-entropy Loss: ' + str(loss_lst.mean())[:7] + ' , ***')
                else:
                    print('Validation at Epoch ' + str(epo + 1) + ', AUROC: ' + str(auc)[:7] + ' , AUPRC: ' + str(auprc)[:7] + ' , F1: ' + str(f1)[:7] + ' , Cross-entropy Loss: ' + str(loss_lst.mean())[:7])

            table.add_row(lst)
        
        #### after training
        prettytable_file = os.path.join(self.result_folder, "valid_markdowntable.txt")
        with open(prettytable_file, 'w') as fp:
            fp.write(table.get_string())
        print('--- Finished ---')

        print('--- Go for Testing ---')
        y_pred = []
        y_label = []
        loss_lst=[]
        model_max.eval()
        for i, (v_d, v_p, label) in enumerate(testing_generator):
            if self.drug_encoding in ["MPNN", 'Transformer', 'DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking','DGL_GIN_ContextPred', 'DGL_AttentiveFP']:
                v_d = v_d
            else:
                v_d = v_d.float().to(self.device)
            if self.target_encoding == 'Transformer':
                v_p = v_p
            else:
                v_p = v_p.float().to(self.device)
            score = model_max(v_d, v_p)
            label = Variable(torch.from_numpy(np.array(label)).float()).to(self.device)
            #m = torch.nn.Sigmoid()
            # ===========loss function1============
            loss_fct = torch.nn.BCELoss()
            m = torch.nn.Sigmoid()
            n = torch.squeeze(m(score), 1)
            loss = loss_fct(n, label)
            # ====================================

            # ===========loss function2============
            #loss_fct = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor(2)) #
            #m = torch.nn.Sigmoid()
            #n = torch.squeeze(score, 1)
            #loss = loss_fct(n, label)
            # ====================================

            loss_lst.append(loss.item())

            logits = torch.squeeze(m(score)).detach().cpu().numpy()
            label_ids = label.to('cpu').numpy()
            y_label = y_label + label_ids.flatten().tolist()
            y_pred = y_pred + logits.flatten().tolist()
            outputs = np.asarray([1 if i else 0 for i in (np.asarray(y_pred) >= 0.5)])
        model_max.train()

        ## ROC-AUC curve
        roc_auc_file = os.path.join(self.result_folder, "roc-auc.jpg")
        plt.figure(0)
        roc_curve(y_pred, y_label, roc_auc_file, self.drug_encoding + '_' + self.target_encoding)
        plt.figure(1)
        pr_auc_file = os.path.join(self.result_folder, "pr-auc.jpg")
        prauc_curve(y_pred, y_label, pr_auc_file, self.drug_encoding + '_' + self.target_encoding)

        ## Output
        auc = roc_auc_score(y_label, y_pred)
        auprc = average_precision_score(y_label, y_pred)
        f1 = f1_score(y_label, outputs)
        #loss = log_loss(y_label, outputs)
        loss_lst = np.array(loss_lst)

        test_table = PrettyTable(["AUROC", "AUPRC", "F1"])
        test_table.add_row(list(map(float2str, [auc, auprc, f1])))

        plt.figure(2)
        cm_file = os.path.join(self.result_folder, "cm.jpg")
        cm = confusion_matrix(y_label, outputs, labels=None, sample_weight=None)
        #print(cm)
        plt.imshow(cm, interpolation='nearest', cmap='Pastel1')
        plt.title('Confusion matrix', size=15)
        plt.colorbar()
        tick_marks = np.arange(2)
        plt.xticks(tick_marks, ["True", "False"], rotation=45, size=10)
        plt.yticks(tick_marks, ["True", "False"] ,size=10)
        plt.tight_layout()
        plt.ylabel('Actual label',size=15)
        plt.xlabel('Predicted label',size=15)
        width, height=cm.shape
        for x in range(width):
            for y in range(height):
                plt.annotate(str(cm[x][y]),xy=(y,x),horizontalalignment='center',verticalalignment='center')
        plt.savefig(cm_file)
        #plt.show()
        #print('Validation at Epoch ' + str(epo + 1) + ' , AUROC: ' + str(auc)[:7] + \' , AUPRC: ' + str(auprc)[:7] + ' , F1: ' + str(f1)[:7] + ' , Cross-entropy Loss: ' + \str(loss)[:7])

        print('Testing at Epoch ' + str(best_epo) + ' , AUROC: ' + str(auc)[:7] + ' , AUPRC: ' + str(auprc)[:7] + ' , F1: ' + str(f1)[:7] + ' , Cross-entropy Loss: ' + str(loss_lst.mean())[:7])

        ######### learning record ###########

        ### 1. test results
        prettytable_file = os.path.join(self.result_folder, "test_markdowntable.txt")
        with open(prettytable_file, 'w') as fp:
            fp.write(test_table.get_string())

        ### 2. learning curve
        fontsize = 16
        iter_num = list(range(1, len(loss_history) + 1))
        plt.figure(3)
        plt.plot(iter_num, loss_history, "bo-")
        plt.xlabel("iteration", fontsize=fontsize)
        plt.ylabel("loss value", fontsize=fontsize)
        pkl_file = os.path.join(self.result_folder, "loss_curve_iter.pkl")
        with open(pkl_file, 'wb') as pck:
            pickle.dump(loss_history, pck)

        fig_file = os.path.join(self.result_folder, "loss_curve.png")
        plt.savefig(fig_file)

        print('--- Training Finished ---')
        writer.flush()
        writer.close()

    """
    def predict(self, df_data):
        '''
            utils.data_process_repurpose_virtual_screening
            pd.DataFrame
        '''
        print('predicting...')
        info = data_process_loader(df_data.index.values, df_data.Label.values, df_data, **self.config)
        self.model.to(device)
        params = {'batch_size': self.config['batch_size'],
                  'shuffle': False,
                  'num_workers': self.config['num_workers'],
                  'drop_last': False,
                  'sampler': SequentialSampler(info)}

        if (self.drug_encoding == "MPNN"):
            params['collate_fn'] = mpnn_collate_func

        generator = data.DataLoader(info, **params)
        y_pred = []
        y_label = []
        self.model.eval()
        for i, (v_d, v_p, label) in enumerate(generator):
            if self.drug_encoding == "MPNN" or self.drug_encoding == 'Transformer':
                v_d = v_d
            else:
                v_d = v_d.float().to(self.device)
            if self.target_encoding == 'Transformer':
                v_p = v_p
            else:
                v_p = v_p.float().to(self.device)
            score = self.model(v_d, v_p)

            m = torch.nn.Sigmoid()
            logits = torch.squeeze(m(score)).detach().cpu().numpy()

            label_ids = label.to('cpu').numpy()
            y_label = y_label + label_ids.flatten().tolist()
            y_pred = y_pred + logits.flatten().tolist()
            result = np.asarray(['Yes' if i else 'No' for i in (np.asarray(y_pred) >= 0.5)])
        self.model.train()
        return y_pred, result
    """


    def predict(self, df_data):
        '''
            utils.data_process_repurpose_virtual_screening
            pd.DataFrame
        '''
        print('predicting...')
        info = data_process_loader(df_data.index.values, df_data.Label.values, df_data, **self.config)
        self.model.to(device)
        params = {'batch_size': self.config['batch_size'],
                  'shuffle': False,
                  'num_workers': self.config['num_workers'],
                  'drop_last': False,
                  'sampler': SequentialSampler(info)}

        if (self.drug_encoding == "MPNN"):
            params['collate_fn'] = mpnn_collate_func
        elif self.drug_encoding in ['DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking', 'DGL_GIN_ContextPred','AttentiveFP']:
            params['collate_fn'] = dgl_collate_func

        generator = data.DataLoader(info, **params)

        y_pred = []
        y_label = []
        self.model.eval()
        for i, (v_d, v_p, label) in enumerate(generator):
            if self.drug_encoding in ["MPNN", 'Transformer', 'DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking','DGL_GIN_ContextPred', 'DGL_AttentiveFP']:
                v_d = v_d
            else:
                v_d = v_d.float().to(self.device)
            if self.target_encoding in ["MPNN", 'Transformer', 'DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking','DGL_GIN_ContextPred', 'DGL_AttentiveFP']:
                v_p = v_p
            else:
                v_p = v_p.float().to(self.device)
            score = self.model(v_d, v_p)

            m = torch.nn.Sigmoid()
            logits = torch.squeeze(m(score)).detach().cpu().numpy()

            label_ids = label.to('cpu').numpy()
            y_label = y_label + label_ids.flatten().tolist()
            y_pred = y_pred + logits.flatten().tolist()
            outputs = np.asarray([1 if i else 0 for i in (np.asarray(y_pred) >= 0.5)])
        self.model.train()
        return y_pred

    def save_model(self, path_dir):
        if not os.path.exists(path_dir):
            os.makedirs(path_dir)
        torch.save(self.model.state_dict(), path_dir + '/model.pt')
        save_dict(path_dir, self.config)

    def load_pretrained(self, path):
        if not os.path.exists(path):
            os.makedirs(path)

        if self.device == 'cuda':
            state_dict = torch.load(path)
        else:
            state_dict = torch.load(path, map_location=torch.device('cpu'))
        # to support training from multi-gpus data-parallel:

        """
        if next(iter(state_dict))[:7] == 'module.':
            # the pretrained model is from data-parallel module
            from collections import OrderedDict
            new_state_dict = OrderedDict()
            for k, v in state_dict.items():
                name = k[7:]  # remove `module.`
                new_state_dict[name] = v
            state_dict = new_state_dict
        """
        self.model.load_state_dict(state_dict)

def model_pretrained(path_dir = None, model = None):
    if model is not None:
        path_dir = download_pretrained_model(model)
    config = load_dict(path_dir)
    model = SmartNet(**config)
    model.load_pretrained(path_dir + '/model.pt')
    return model

def classify(X_repurpose, target, model, drug_names=None, target_names=None, result_folder="smartnet/smartnet_vscode/result/", output_num_max=10):
    # X_repurpose: a list of SMILES string
    # target: a list of targets 


    fo = os.path.join(result_folder, "Results.txt")
    print_list = []

    if drug_names is None:
        drug_names = ['Drug ' + str(i) for i in list(range(len(X_repurpose)))]
    if target_names is None:
        target_names = ['Target ' + str(i) for i in list(range(len(target)))]

    table_header = ["Rank", "Drug Name", "Target Name", "Interaction", "Probability"]
    table = PrettyTable(table_header)



    with open(fo, 'w') as fout:
        print('virtual screening...')
        df_data = data_process_repurpose_virtual_screening(X_repurpose, target, model.drug_encoding, model.target_encoding)

        #print(df_data)
        y_pred = model.predict(df_data)
        #print(y_pred)

        print('---------------')
        """
        if drug_names is not None and target_names is not None:
            print('Virtual Screening Result')
            f_d = max([len(o) for o in drug_names]) + 1
            f_p = max([len(o) for o in target_names]) + 1
        """

        #print(df_data)
        #print(len(X_repurpose), len(drug_names), len(target_names))

        #print(len(target_names))
        #if isinstance(target_names, str)==True:
        d_n = []
        t_n = []
        if len(target_names) == 1:
            t_n = np.tile(target_names, (len(drug_names),1)).reshape(-1,).tolist()
            d_n=drug_names

        elif len(drug_names) == 1:
            d_n = np.tile(drug_names, (len(target_names),1)).reshape(-1,).tolist()
            t_n=target_names
        elif len(target_names) >1 and len(drug_names) > 1:
            for i in range(len(drug_names)):
                for j in range(len(target_names)):
                    d_n.append(drug_names[i])
                    t_n.append(target_names[j])


        #print(d_n)
        #print(t_n)
        #print(len(X_repurpose), len(drug_names), len(target_names))
        #print(type(target_names))
        #print(target_names)

        #print(len(target_names))
        #print(len(t_n))
        #print(len(d_n))


        for i in range(len(d_n)):
            if y_pred[i] > 0.5:
                string_lst = [d_n[i], t_n[i], "YES", "{0:.2f}".format(y_pred[i])]

            else:
                string_lst = [d_n[i], t_n[i], "NO", "{0:.2f}".format(y_pred[i])]

            print_list.append((string_lst, y_pred[i]))


        print_list.sort(key=lambda x: x[1], reverse=True)

        print_list = [i[0] for i in print_list]
        for idx, lst in enumerate(print_list):
            lst = [str(idx + 1)] + lst
            table.add_row(lst)
        fout.write(table.get_string())


    with open(fo, 'r') as fin:
        lines = fin.readlines()
        for idx, line in enumerate(lines):
            if idx < output_num_max+3:
                print(line, end='')
            else:
                print('checkout ' + fo + ' for the whole list')
                break

    return y_pred

def classify_ss(X_repurpose, target, structure, model, drug_names=None, target_names=None, result_folder="smartnet/smartnet_vscode/result", output_num_max=10):
    # X_repurpose: a list of SMILES string
    # target: a list of targets


    fo = os.path.join(result_folder, "Results.txt")
    print_list = []

    if drug_names is None:
        drug_names = ['Drug ' + str(i) for i in list(range(len(X_repurpose)))]
    if target_names is None:
        target_names = ['Target ' + str(i) for i in list(range(len(target)))]

    table_header = ["Rank", "Drug Name", "Target Name", "Interaction", "Probability"]
    table = PrettyTable(table_header)



    with open(fo, 'w') as fout:
        print('virtual screening...')
        df_data = data_process_repurpose_virtual_screening_ss(X_repurpose, target, structure, model.drug_encoding, model.target_encoding)

        #print(df_data)
        y_pred = model.predict(df_data)
        #print(y_pred)

        print('---------------')
        """
        if drug_names is not None and target_names is not None:
            print('Virtual Screening Result')
            f_d = max([len(o) for o in drug_names]) + 1
            f_p = max([len(o) for o in target_names]) + 1
        """

        #print(df_data)
        #print(len(X_repurpose), len(drug_names), len(target_names))

        #print(len(target_names))
        #if isinstance(target_names, str)==True:
        d_n = []
        t_n = []
        if len(target_names) == 1:
            t_n = np.tile(target_names, (len(drug_names),1)).reshape(-1,).tolist()
            d_n=drug_names

        elif len(drug_names) == 1:
            d_n = np.tile(drug_names, (len(target_names),1)).reshape(-1,).tolist()
            t_n=target_names
        elif len(target_names) >1 and len(drug_names) > 1:
            for i in range(len(drug_names)):
                for j in range(len(target_names)):
                    d_n.append(drug_names[i])
                    t_n.append(target_names[j])


        #print(d_n)
        #print(t_n)
        #print(len(X_repurpose), len(drug_names), len(target_names))
        #print(type(target_names))
        #print(target_names)

        #print(len(target_names))
        #print(len(t_n))
        #print(len(d_n))


        for i in range(len(d_n)):
            if y_pred[i] > 0.5:
                string_lst = [d_n[i], t_n[i], "YES", "{0:.2f}".format(y_pred[i])]

            else:
                string_lst = [d_n[i], t_n[i], "NO", "{0:.2f}".format(y_pred[i])]

            print_list.append((string_lst, y_pred[i]))


        print_list.sort(key=lambda x: x[1], reverse=True)

        print_list = [i[0] for i in print_list]
        for idx, lst in enumerate(print_list):
            lst = [str(idx + 1)] + lst
            table.add_row(lst)
        fout.write(table.get_string())


    with open(fo, 'r') as fin:
        lines = fin.readlines()
        for idx, line in enumerate(lines):
            if idx < output_num_max+3:
                print(line, end='')
            else:
                print('checkout ' + fo + ' for the whole list')
                break

    return y_pred