#!/usr/bin/env python
# coding: utf-8

# In[1]:


#%%
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import loader
import processor
import config
import controller
import models


# In[2]:


from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))


# In[79]:


loader.valid_ss(input='./data/final_data1.txt',output='./data/final_data1_valid.txt')
#loader.valid(input='./data/final_data0.txt',output='./data/final_data0_valid.txt')


# # 固定随机种子

# In[80]:


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


# In[81]:


fix_seed(1)


# # 数据分割方式改变

# In[82]:


drug_encoding, target_encoding = 'CNN', 'CNN'


# In[83]:


train, val, test=loader.split_drug(input='./data/final_data1_valid.txt', random_seed=1)


# In[84]:


train


# # 根据1:2的比例等比例划分数据集

# In[85]:


X_drugs_train, X_targets_train, X_structures_train, y_train=loader.df2var_bal_ss(train,weight=2, random_seed=10000)
#X_drugs_train, X_targets_train, y_train=loader.df2var_bal(train,weight=2, random_seed=10000)


# In[86]:


X_drugs_val, X_targets_val,  X_structures_val, y_val=loader.df2var_bal_ss(val,weight=2, random_seed=20000)
#X_drugs_val, X_targets_val, y_val=loader.df2var_bal(val,weight=2, random_seed=20000)


# In[87]:


X_drugs_test, X_targets_test, X_structures_test, y_test=loader.df2var_bal_ss(test,weight=2, random_seed=30000)
#X_drugs_test, X_targets_test, y_test=loader.df2var_bal(test,weight=2, random_seed=30000)


# In[ ]:





# In[ ]:





# In[88]:


train = processor.encode_ss(X_drugs_train, X_targets_train, X_structures_train, y_train, drug_encoding, target_encoding, random_seed = 1)
#train = processor.encode(X_drugs_train, X_targets_train, y_train, drug_encoding, target_encoding, random_seed = 1)


# In[89]:


train


# In[ ]:





# In[90]:


val = processor.encode_ss(X_drugs_val, X_targets_val, X_structures_val, y_val, drug_encoding, target_encoding, random_seed = 1)
#val = processor.encode(X_drugs_val, X_targets_val, y_val, drug_encoding, target_encoding, random_seed = 1)


# In[91]:


test = processor.encode_ss(X_drugs_test, X_targets_test, X_structures_test, y_test, drug_encoding, target_encoding, random_seed = 1)
#test = processor.encode(X_drugs_test, X_targets_test, y_test, drug_encoding, target_encoding, random_seed = 1)


# In[22]:


#train = train.loc[:,train.columns!='structure_encoding']
#train = train.loc[:,train.columns!='Target Structure']
#val = val.loc[:,val.columns!='structure_encoding']
#val = val.loc[:,val.columns!='Target Structure']
#test = test.loc[:,test.columns!='structure_encoding']
#test = test.loc[:,test.columns!='Target Structure']


# In[92]:


train


# In[93]:


val


# In[94]:


test


# In[95]:


import torch.nn.functional as F

param = config.set(drug_encoding, 
                   target_encoding, 
                   result_folder = "./result/",
                   #input_dim_drug = 1024, 
                   #input_dim_protein = 8420,
                   hidden_dim_drug = 128, 
                   hidden_dim_protein = 128,
                   cls_hidden_dims = [1024,1024,512], 
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


# In[96]:


model = controller.initialize(**param)


# In[97]:


model


# In[98]:


#model_test = models.Prism_GIN_AttrMasking(**param)
model_test = models.CNN_CNN(**param)


# In[99]:


print(model_test)


# In[100]:


#text_fixed #20210819 no shffule
train


# In[101]:


model.train(train, val, test)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


#text_fixed #20210819 shuffle


# In[22]:


model.train(train, val, test)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[24]:


#text_fixed2


# In[27]:


model.train(train, val, test)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[67]:


model.train(train, val, test)


# In[ ]:





# In[ ]:


#retest


# In[19]:


model.train(train, val, test)


# In[ ]:





# In[ ]:


#


# In[19]:


model.train(train, val, test)


# In[ ]:





# In[19]:


#20200805 retest


# In[19]:


model.train(train, val, test)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


seed=1, free, 4298/533/835


# In[17]:


model.train(train, val, test)


# In[ ]:





# In[ ]:


seed=1, seed2=1,1,1 


# In[19]:


model.train(train, val, test)


# In[ ]:





# In[ ]:





# In[ ]:


seed1=1, seed2-10000,20000,30000


# In[18]:


model.train(train, val, test)


# In[ ]:


6:2:2


# In[18]:


model.train(train, val, test)


# In[ ]:


seed1=1, seed2-10000,20000,30000, 7:2:1


# In[18]:


model.train(train, val, test)


# In[ ]:





# In[ ]:


seed1=1, seed2-10000,20000,30000, 6:2:2


# In[35]:


model.train(train, val, test)


# In[ ]:





# In[ ]:


seed1=1, seed2-10000,20000,30000, 8:1:1


# In[18]:


model.train(train, val, test)


# In[ ]:





# In[ ]:





# In[20]:


model.save_model('./Prism_GIN_drug_bal_model_seq')


# In[22]:


model = controller.model_pretrained(path_dir = './Prism_GIN_drug_bal_model_seq')
model


# In[49]:


t_name, t= loader.predict_load('./case_study/ESE2_target.txt', 1)


# In[50]:


d_name, d= loader.predict_load('./case_study/ESE2_drug.txt', 1)


# In[51]:


y_pred = controller.classify(X_repurpose = d, 
                             target = t, 
                             model = model, 
                             drug_names = d_name, 
                             target_names = t_name, 
                             result_folder = "./result/", 
                             output_num_max=80)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[58]:


t_name, t= loader.predict_load('./case_study/disney_target.txt', 1)


# In[59]:


d_name, d= loader.predict_load('./case_study/disney_drug.txt', 5)


# In[60]:


y_pred = controller.classify(X_repurpose = d, 
                             target = t, 
                             model = model, 
                             drug_names = d_name, 
                             target_names = t_name, 
                             result_folder = "./result/", 
                             output_num_max=80)


# In[ ]:





# In[ ]:





# In[ ]:





# In[17]:


model


# In[21]:


#%%
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import os
import datautils
from smoothgrad import *
from param import param_num


# In[147]:


X_drugs, X_targets, y = loader.file2var(input = './data/gra_test.txt')


# In[148]:


drug_encoding, target_encoding = 'DGL_GIN_AttrMasking', 'CNN'


# In[149]:


test = processor.encode(X_drugs, X_targets, y, drug_encoding, target_encoding, random_seed = 1)


# In[150]:


import torch.nn.functional as F

param = config.set(drug_encoding, 
                   target_encoding, 
                   result_folder = "./result/",
                   #input_dim_drug = 1024, 
                   #input_dim_protein = 8420,
                   hidden_dim_drug = 128, 
                   hidden_dim_protein = 128,
                   cls_hidden_dims = [1024,1024,512], 
                   #batch_size = 256, 
                   batch_size = 64, 
                   train_epoch = 5, 
                   LR = 0.001, 
                   cnn_drug_filters = [32,64,96],
                   cnn_drug_kernels = [4,6,8], #odd
                   cnn_target_filters = [32,64,96],
                   cnn_target_kernels = [4,6,8], #odd
                   gnn_hid_dim_drug = 64,
                   gnn_num_layers = 3,
                   gnn_activation = F.relu,
                   in_feats = 74
                  )


# In[151]:


import dgl
param['batch_size']=1
def dgl_collate_func(x):
    d, p, y = zip(*x)
    d = dgl.batch(d)
    return d, torch.tensor(p), torch.tensor(y)

import torch
from torch.utils.data import SequentialSampler
from loader import *
params_test = {'batch_size': param['batch_size'],
               'shuffle': False,
               'num_workers': 0,
               'drop_last': False,
                'collate_fn': dgl_collate_func,
               'sampler': SequentialSampler(data_process_loader(test.index.values, test.Label.values, test, **param))}

test_loader = torch.utils.data.DataLoader(data_process_loader(test.index.values, test.Label.values, test, **param), **params_test)


# In[152]:


import torch
device = torch.device("cuda")


# In[153]:


model = models.CNN_GIN_AttrMasking(**param)


# In[154]:


model_path = ("./CNN_GIN_drug_bal_model/model.pt")
filename = model_path.format("CNN_GIN")
print("Loading model: {}".format(filename))
model.load_state_dict(torch.load(filename, map_location='cpu'))
model = model.to(device)
#print(type(model))


# In[155]:


from smoothgrad import *

batch_size=1

def compute_saliency(model, device, test_loader):

    model.eval()

    identity = "result"
    #saliency_dir = datautils.make_directory(".", "out/saliency")
    saliency_path = os.path.join("./out", identity+'.sal')

    # sgrad = SmoothGrad(model, device=device)
    sgrad = GuidedBackpropSmoothGrad(model, device=device)
    sal = ""
    #for batch_idx, (x0, y0) in enumerate(test_loader):
        #X, Y = x0.float().to(device), y0.to(device).float()
    #print(len(test))
    for batch_idx, (v_d, v_p, label) in enumerate(test_loader):
        X, Y, Z = v_d, v_p.float().to(device), label.to(device)
        #print(X.shape, Y.shape)
        #print(batch_idx)
        guided_saliency_p= sgrad.get_batch_gradients(X, Y, Z)
        # import pdb; pdb.set_trace()
        N, NS, _ = guided_saliency_p.shape # (N, 31, 1, 4)
        
        output = model(X, Y)
        prob = torch.sigmoid(output)
        p_np = prob.to(device='cpu').detach().numpy().squeeze()
        
        #print(N)
        str_sal=[] 
        for i in range(N):
            inr = batch_idx*batch_size + i
            str_sal = datautils.mat2str(np.squeeze(guided_saliency_p[i]))
            #print(p_np)
            #sal += "{}\t{:.6f}\t{}\n".format(inr, p_np, str_sal)
            str_sal
    return str_sal[:-1]
       
    #f = open(saliency_path,"w")
    #f.write(sal)
    #f.close()
    #print(saliency_path)


# In[156]:


sal=compute_saliency(model, device, test_loader)


# In[157]:


sal


# In[158]:


sal_list=sal.split(",")


# In[159]:


test['Target Sequence'][0]


# In[160]:


seq_list=list(test['Target Sequence'][0])


# In[161]:


seq = ','.join(seq_list)


# In[162]:


print(seq)
for n in range(len(sal_list)//31):
    print(','.join(sal_list[n*31:(n+1)*31]))


# In[163]:


with open('./data/df.txt', 'w+') as f:
    #f.write(str(seq) + '\n')
    for n in range(len(sal_list)//31):
        f.write(','.join(sal_list[n*31:(n+1)*31]) + '\n')


# In[164]:


df=pd.read_table('./data/df.txt',sep=',',header=None)


# In[165]:


df


# In[166]:


df2=pd.read_table('./data/df2.txt',sep=',',header=None)


# In[167]:


import pandas as pd
df3 = pd.Series(list('AAAGGUUUAUACCUUCCCAGGUAACAAACCA'))
df4=(pd.get_dummies(df3)).T


# In[168]:


df5=df4.rename(index={'A': '0','U': '1','C': '2','G': '3'})


# In[169]:


df6=df5.sort_index()


# In[170]:


merge = df6.append(df2, ignore_index=True)


# In[171]:


df_final=df*merge


# In[172]:


df_final


# In[173]:


df_final.max().max()


# In[174]:


df_norm=df_final/df_final.max().max()


# In[175]:


import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns


# In[176]:


plt.figure(figsize=(30, 5))

h=sns.heatmap(df_norm, cmap='Reds', linewidths=0.1, linecolor='grey',annot=False,cbar=False)
cb = h.figure.colorbar(h.collections[0]) #显示colorbar
cb.ax.tick_params(labelsize=10)  # 设置colorbar刻度字体大小。
#,cbar_kws={"shrink": 0.5}
#plt.savefig("chr.pdf")
plt.show()


# In[ ]:





# In[71]:


df_all=pd.read_table('./data/df.txt',sep=',',header=None)


# In[72]:


df_all.max().max()


# In[73]:


df_all_norm=df_all/df_all.max().max()


# In[74]:


df_all_norm


# In[75]:


plt.figure(figsize=(30, 5))

h=sns.heatmap(df_all_norm, cmap='Reds', linewidths=0.1, linecolor='grey',annot=False,cbar=False)
cb = h.figure.colorbar(h.collections[0]) #显示colorbar
cb.ax.tick_params(labelsize=10)  # 设置colorbar刻度字体大小。
#,cbar_kws={"shrink": 0.5}
#plt.savefig("chr.pdf")
plt.show()


# In[ ]:





# In[81]:


t_name, t= loader.predict_load('./data/target2.txt', 1)


# In[82]:


d_name, d= loader.predict_load('./data/drug.txt', 80)


# In[83]:


model = controller.model_pretrained(path_dir = './CNN_GIN_drug_bal_model')
model


# In[84]:


y_pred = controller.classify(X_repurpose = d, 
                             target = t, 
                             model = model, 
                             drug_names = d_name, 
                             target_names = t_name, 
                             result_folder = "./result/", 
                             output_num_max=80)


# In[ ]:


#drug visualization


# In[6]:


from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG


# In[16]:


smiles='CC(=O)OC1=CC=CC=C1C(=O)NC2=NC=C(S2)[N+](=O)[O-]'


# In[17]:


m = Chem.MolFromSmiles(smiles)


# In[18]:


m


# In[19]:


def moltosvg(mol, molSize = (300,300), kekulize = True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:','')

SVG(moltosvg(m))


# In[ ]:


#disney test


# In[7]:


model = controller.model_pretrained(path_dir = './CNN_GIN_drug_bal_model')
model


# In[14]:


t_name, t= loader.predict_load('./data/mir21_target.txt', 3)


# In[15]:


d_name, d= loader.predict_load('./data/mir21_drug.txt', 3)


# In[16]:


y_pred = controller.classify(X_repurpose = d, 
                             target = t, 
                             model = model, 
                             drug_names = d_name, 
                             target_names = t_name, 
                             result_folder = "./result/", 
                             output_num_max=10)


# In[75]:


C1='COc1ccc(Br)cc1C(=O)Nc4ccc(c3nc2ccccc2[nH]3)cc4'


# In[ ]:


C2


# In[ ]:


C3


# In[ ]:


C4


# In[ ]:


C5


# In[83]:


C5='CNc2nc(Nc1cccc(C(=O)O)c1)nc3ccccc23'


# In[47]:


sm1='COC(=O)CCCOc5ccc(c4nc3ccc(c2ccc(N1CCN(C)CC1)cc2)cc3[nH]4)cc5'


# In[48]:


from rdkit import Chem
Chem.MolFromSmiles(C5)


# In[76]:


Chem.MolFromSmiles(C1)


# In[ ]:





# In[ ]:





# In[20]:


#%%
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import os
import datautils
from smoothgrad import *
from param import param_num


# In[21]:


X_drugs, X_targets, y = loader.file2var(input = './data/gra_test.txt')


# In[22]:


drug_encoding, target_encoding = 'DGL_GIN_AttrMasking', 'CNN'


# In[23]:


test = processor.encode(X_drugs, X_targets, y, drug_encoding, target_encoding, random_seed = 1)


# In[24]:


import torch.nn.functional as F

param = config.set(drug_encoding, 
                   target_encoding, 
                   result_folder = "./result/",
                   #input_dim_drug = 1024, 
                   #input_dim_protein = 8420,
                   hidden_dim_drug = 128, 
                   hidden_dim_protein = 128,
                   cls_hidden_dims = [1024,1024,512], 
                   #batch_size = 256, 
                   batch_size = 64, 
                   train_epoch = 5, 
                   LR = 0.001, 
                   cnn_drug_filters = [32,64,96],
                   cnn_drug_kernels = [4,6,8], #odd
                   cnn_target_filters = [32,64,96],
                   cnn_target_kernels = [4,6,8], #odd
                   gnn_hid_dim_drug = 64,
                   gnn_num_layers = 3,
                   gnn_activation = F.relu,
                   in_feats = 74
                  )


# In[25]:


import dgl
param['batch_size']=1
def dgl_collate_func(x):
    d, p, y = zip(*x)
    d = dgl.batch(d)
    return d, torch.tensor(p), torch.tensor(y)

import torch
from torch.utils.data import SequentialSampler
from loader import *
params_test = {'batch_size': param['batch_size'],
               'shuffle': False,
               'num_workers': 0,
               'drop_last': False,
                'collate_fn': dgl_collate_func,
               'sampler': SequentialSampler(data_process_loader(test.index.values, test.Label.values, test, **param))}

test_loader = torch.utils.data.DataLoader(data_process_loader(test.index.values, test.Label.values, test, **param), **params_test)


# In[26]:


import torch
device = torch.device("cuda")


# In[27]:


model = models.CNN_GIN_AttrMasking(**param)


# In[28]:


model_path = ("./CNN_GIN_drug_bal_model/model.pt")
filename = model_path.format("CNN_GIN")
print("Loading model: {}".format(filename))
model.load_state_dict(torch.load(filename, map_location='cuda'))
model = model.to(device)
#print(type(model))


# In[29]:


from smoothgrad import *

batch_size=1

def compute_saliency(model, device, test_loader):

    model.eval()

    identity = "result"
    #saliency_dir = datautils.make_directory(".", "out/saliency")
    saliency_path = os.path.join("./out", identity+'.sal')

    # sgrad = SmoothGrad(model, device=device)
    sgrad = GuidedBackpropSmoothGrad(model, device=device)
    sal = ""
    #for batch_idx, (x0, y0) in enumerate(test_loader):
        #X, Y = x0.float().to(device), y0.to(device).float()
    #print(len(test))
    for batch_idx, (v_d, v_p, label) in enumerate(test_loader):
        X, Y, Z = v_d, v_p.float().to(device), label.to(device)
        #print(X.shape, Y.shape)
        #print(batch_idx)
        guided_saliency_p= sgrad.get_batch_gradients(X, Y, Z)
        # import pdb; pdb.set_trace()
        N, NS, _ = guided_saliency_p.shape # (N, 31, 1, 4)
        
        output = model(X, Y)
        prob = torch.sigmoid(output)
        p_np = prob.to(device='cpu').detach().numpy().squeeze()
        
        #print(N)
        str_sal=[] 
        for i in range(N):
            inr = batch_idx*batch_size + i
            str_sal = datautils.mat2str(np.squeeze(guided_saliency_p[i]))
            #print(p_np)
            #sal += "{}\t{:.6f}\t{}\n".format(inr, p_np, str_sal)
            str_sal
    return str_sal[:-1]
       
    #f = open(saliency_path,"w")
    #f.write(sal)
    #f.close()
    #print(saliency_path)


# In[30]:


sal=compute_saliency(model, device, test_loader)


# In[31]:


sal_list=sal.split(",")


# In[32]:


test['Target Sequence'][0]
seq_list=list(test['Target Sequence'][0])
seq = ','.join(seq_list)
with open('./data/df.txt', 'w+') as f:
    #f.write(str(seq) + '\n')
    for n in range(len(sal_list)//31):
        f.write(','.join(sal_list[n*31:(n+1)*31]) + '\n')
df=pd.read_table('./data/df.txt',sep=',',header=None)
df2=pd.read_table('./data/df2.txt',sep=',',header=None)
import pandas as pd
df3 = pd.Series(list('AAUCUCAUGGCAACACCAGUCGAUGGGCUGU'))
df4=(pd.get_dummies(df3)).T
df5=df4.rename(index={'A': '0','U': '1','C': '2','G': '3'})
df6=df5.sort_index()
merge = df6.append(df2, ignore_index=True)
df_final=df*merge
df_final
df_final.max().max()
df_norm=df_final/df_final.max().max()
import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns
plt.figure(figsize=(30, 5))
h=sns.heatmap(df_norm, cmap='Reds', linewidths=0.1, linecolor='grey',annot=False,cbar=False)
cb = h.figure.colorbar(h.collections[0]) #显示colorbar
cb.ax.tick_params(labelsize=10)  # 设置colorbar刻度字体大小。
#,cbar_kws={"shrink": 0.5}
#plt.savefig("chr.pdf")
plt.show()


# In[ ]:




