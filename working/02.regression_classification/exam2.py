import loader
import processor
import config
import controller
import models
import utils

#inputpath1 = 'smartnet/smartnet_vscode/data/final_data_autodock_knot1.txt'
#outputpath1 = 'smartnet/smartnet_vscode/data/final_data_autodock_knot_valid.txt'
#outputpath2 = 'smartnet/smartnet_vscode/result/regression/GIN_Prism_autodock_SPUK2'
#outputpath3 = 'smartnet/smartnet_vscode/result/regression/GIN_Prism_autodock_SPUK2/GIN_Prism_autodock_SPUK2'

loader.valid_ss(input='smartnet/smartnet_classification/data/final_data1.txt',output='smartnet/smartnet_classification/data/final_data1_valid.txt')
#loader.valid(input='smartnet/smartnet_classification/data/final_data0.txt',output='smartnet/smartnet_classification/data/final_data0_valid.txt')

import random
import torch
import os
import numpy as np
def fix_seed(seed):
    """
    Seed all necessary random number generators.
    """
    #if seed is None:
    #    seed = random.randint(1, 10000)
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

fix_seed(1)

drug_encoding, target_encoding = 'CNN', 'CNN'
#drug_encoding, target_encoding = 'DGL_GIN_AttrMasking','Prism'
#drug_encoding, target_encoding = 'DGL_GIN_AttrMasking','CNN'
#X_drugs, X_targets, X_structures,X_knots, y =loader.file2var_reg_ss2(inputpath1, convert_to_log = False)
#X_drugs, X_targets, X_structures, y =loader.file2var_reg_ss(inputpath1, convert_to_log = False)
train, val, test=loader.split_drug(input='smartnet/smartnet_classification/data/final_data1_valid.txt', random_seed=1)

# # 根据1:2的比例等比例划分数据集
"""
train, val, test = utils.data_process_ss2(X_drugs, X_targets, X_structures,X_knots,
                                y, drug_encoding, target_encoding, 
                                split_method='random',frac=[0.7,0.1,0.2],
                                random_seed = 1)
"""
X_drugs_train, X_targets_train, X_structures_train, y_train=loader.df2var_bal_ss(train,weight=2, random_seed=10000)
#X_drugs_train, X_targets_train, y_train=loader.df2var_bal(train,weight=2, random_seed=10000)

X_drugs_val, X_targets_val,  X_structures_val, y_val=loader.df2var_bal_ss(val,weight=2, random_seed=20000)
#X_drugs_val, X_targets_val, y_val=loader.df2var_bal(val,weight=2, random_seed=20000)

X_drugs_test, X_targets_test, X_structures_test, y_test=loader.df2var_bal_ss(test,weight=2, random_seed=30000)
#X_drugs_test, X_targets_test, y_test=loader.df2var_bal(test,weight=2, random_seed=30000)

train = processor.encode_ss(X_drugs_train, X_targets_train, X_structures_train, y_train, drug_encoding, target_encoding, random_seed = 1)
#train = processor.encode(X_drugs_train, X_targets_train, y_train, drug_encoding, target_encoding, random_seed = 1)

val = processor.encode_ss(X_drugs_val, X_targets_val, X_structures_val, y_val, drug_encoding, target_encoding, random_seed = 1)
#val = processor.encode(X_drugs_val, X_targets_val, y_val, drug_encoding, target_encoding, random_seed = 1)

test = processor.encode_ss(X_drugs_test, X_targets_test, X_structures_test, y_test, drug_encoding, target_encoding, random_seed = 1)
#test = processor.encode(X_drugs_test, X_targets_test, y_test, drug_encoding, target_encoding, random_seed = 1)

#train = train.loc[:,train.columns!='structure_encoding']
#train = train.loc[:,train.columns!='Target Structure']
#val = val.loc[:,val.columns!='structure_encoding']
#val = val.loc[:,val.columns!='Target Structure']
#test = test.loc[:,test.columns!='structure_encoding']
#test = test.loc[:,test.columns!='Target Structure']
#train
#val
#test

import torch.nn.functional as F
param = config.set(drug_encoding, 
                   target_encoding, 
                   result_folder = 'smartnet/smartnet_vscode/result/CNN_CNN_inner', #outputpath2,
                   #input_dim_drug = 1024, 
                   #input_dim_protein = 8420,
                   hidden_dim_drug = 128, 
                   hidden_dim_protein = 128,
                   cls_hidden_dims = [1024,512,256], 
                   #batch_size = 256, 
                   batch_size = 64, 
                   train_epoch = 100, 
                   LR = 0.001, 
                   cnn_drug_filters = [32,64,96],
                   cnn_drug_kernels = [4,6,8], #odd
                   cnn_target_filters = [32,64,96],
                   cnn_target_kernels = [4,6,8], #odd
                   num_workers=0
                  )

model = controller.initialize(**param)
#model_test = models.CNN_CNN_mul(**param)
#model_test = models.CNN_GIN_AttrMasking_mul(**param)
#model_test = models.Prism_GIN_AttrMasking_mul(**param)
#model_test = models.Prism_GIN_AttrMasking_TF(**param)
model_test = models.CNN_CNN_mul_DTI(**param)
print(model_test)
print(train)

model.train(train, val, test)
model.save_model('smartnet/smartnet_vscode/result/CNN_CNN_inner/CNN_CNN_drug_bal_model_seq')
#model.train(train, val, test)
