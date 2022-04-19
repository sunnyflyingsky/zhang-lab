import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
 
from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA


def prepareData(flag='drug'):
    '''
    数据预处理，去重以及合并等
    '''
    if flag=='drug':
        ref_FDA = {}
        with open('docking/data/FDA_valid.txt') as read_object:
            for line in read_object:
                info = line.strip().split('\t')
                ref_FDA[info[1]]=info[0]
        print(len(ref_FDA))
        with open('docking/data/FDA_valid_1.txt','w') as write_object:
            for key,value in ref_FDA.items():
                write_object.write('{}\t{}\n'.format(value,key))

        
        ref_valid = {}
        with open('docking/data/PDB_valid.txt') as read_object:
            for line in read_object:
                info = line.strip().split('\t')
                ref_valid[info[1]]=info[0]
        print(len(ref_valid))
        with open('docking/data/PDB_valid_1.txt','w') as write_object:
            for key,value in ref_valid.items():
                write_object.write('{}\t{}\n'.format(value,key))

        with open('docking/data/ALL_valid_1.txt','w') as write_object:
            for key,value in ref_valid.items():
                write_object.write('{}\t{}\t{}\n'.format(value,key,1))
            for key,value in ref_FDA.items():
                if key not in ref_valid.keys():
                    ref_valid[key] = value
                    write_object.write('{}\t{}\t{}\n'.format(value,key,0))
            print(len(ref_valid))

        ref_valid = {}
        ref_filtered = {}
        with open('docking/data/ALL_valid_1.txt') as read_object:
            for line in read_object:
                info = line.strip().split('\t')
                ref_valid[info[1]]=(info[0],int(info[2]))
        write_object1 = open('docking/data/ALL_valid_2.txt','w')
        write_object2 = open('docking/data/PDB_valid_2.txt','w')
        write_object3 = open('docking/data/FDA_valid_2.txt','w')
        with open('docking/data/results_PDB.txt') as read_object:
            for line in read_object:
                smile = line.strip().split('\t')[0]
                if smile not in ref_filtered.keys():
                    write_object1.write('{}\t{}\t{}\n'.format(ref_valid[smile][0],smile,ref_valid[smile][1]))
                    if ref_valid[smile][1]:
                        write_object2.write('{}\t{}\t{}\n'.format(ref_valid[smile][0],smile,ref_valid[smile][1]))
                    else:
                        write_object3.write('{}\t{}\t{}\n'.format(ref_valid[smile][0],smile,ref_valid[smile][1]))
                    ref_filtered[smile]=1
        with open('docking/data/results_docking.txt') as read_object:
            for line in read_object:
                smile = line.strip().split('\t')[0]
                if smile not in ref_filtered.keys():
                    write_object1.write('{}\t{}\t{}\n'.format(ref_valid[smile][0],smile,ref_valid[smile][1]))
                    if ref_valid[smile][1]:
                        write_object2.write('{}\t{}\t{}\n'.format(ref_valid[smile][0],smile,ref_valid[smile][1]))
                    else:
                        write_object3.write('{}\t{}\t{}\n'.format(ref_valid[smile][0],smile,ref_valid[smile][1]))
                    ref_filtered[smile]=1
    elif flag=='RNA':
        ref_filtered = {}
        write_object1 = open('docking/data/ALL_valid_3.txt','w')
        write_object2 = open('docking/data/PDB_valid_3.txt','w')
        write_object3 = open('docking/data/FDA_valid_3.txt','w')
        i = 0
        with open('docking/data/results_PDB.txt') as read_object:
            for line in read_object:
                info = line.strip().split('\t')
                seq = info[1]
                struct = info[2]
                if seq+'@'+struct not in ref_filtered.keys():
                    i+=1
                    write_object1.write('{}\t{}\t{}\t{}\n'.format(str(i),seq,struct,1))
                    write_object2.write('{}\t{}\t{}\t{}\n'.format(str(i),seq,struct,1))
                    ref_filtered[seq+'@'+struct]=1
        with open('docking/data/results_docking.txt') as read_object:
            for line in read_object:
                info = line.strip().split('\t')
                seq = info[1]
                struct = info[2]
                if seq+'@'+struct not in ref_filtered.keys():
                    i+=1
                    write_object1.write('{}\t{}\t{}\t{}\n'.format(str(i),seq,struct,0))
                    write_object3.write('{}\t{}\t{}\t{}\n'.format(str(i),seq,struct,0))
                    ref_filtered[seq+'@'+struct]=1
    else:
        pass

def calcDistribution(name='FDA'):
    '''
    计算碱基分布曲线
    '''
    prob = {'A':[0]*31,'U':[0]*31,'C':[0]*31,'G':[0]*31,'N':[0]*31}
    j = 0
    with open('docking/data/'+name+'_valid_3.txt') as read_object:
        for line in read_object:
            seq = line.strip().split('\t')[1]
            for i in range(len(seq)):
                prob[seq[i:i+1]][i]+=1
            j+=1
    for key,value in prob.items():
        for i in range(len(value)):
            prob[key][i] /= j
    plt.figure(figsize=(10, 5),dpi=300)
    plt.plot(list(range(1,32,1)),prob['A'],c='red',label='A')
    plt.plot(list(range(1,32,1)),prob['U'],c='blue',label='U')
    plt.plot(list(range(1,32,1)),prob['C'],c='yellow',label='C')
    plt.plot(list(range(1,32,1)),prob['G'],c='green',label='G')
    plt.plot(list(range(1,32,1)),prob['N'],c='black',label='N')
    plt.title('Nucleotide Distribution')
    plt.legend()
    plt.savefig('docking/cluster_res/distrib/'+name+'.png')
            
def loadingDrug(inputpath):
    '''
    加载药物
    '''
    data = []
    with open(inputpath) as read_object:
        i = 1
        for line in read_object:
            info = line.strip().split('\t')
            molecule = info[0]
            smile = info[1]
            label = int(info[2])
            #molecule = 'compound'+str(i)
            #smile = info[0]
            mol = Chem.MolFromSmiles(smile)    # 转换SMILES 为rdkit分子对象 
            mol.SetProp('_Name',molecule)       # 为每个分子添加名字
            #fps = FingerprintMols.FingerprintMol(mol)
            fps = Chem.MACCSkeys.GenMACCSKeys(mol)
            data.append([molecule,smile,mol,fps,label])
            i+=1
    print('loading '+str(i)+' drugs over!')
    return data

def loadingRNA(inputpath):
    '''
    加载RNA
    '''
    database = []
    with open(inputpath) as read_object:
        for line in read_object:
            info = line.strip().split('\t')
            database.append(info)
    return database

def calcSmilarity(fps1,fps2):
    '''
    计算两个化合物小分子的相似性
    '''
    return DataStructs.FingerprintSimilarity(fps1,fps2)

def calcMatrix(database):
    '''
    计算小分子两两之间的相似性
    '''
    size=len(database)
    hmap=np.empty(shape=(size,size))
    table=pd.DataFrame()
    for index, molObjecti in enumerate(database):
        print(index)
        for jndex, molObjectj in enumerate(database):
            similarity=calcSmilarity(molObjecti[3],molObjectj[3])
            hmap[index,jndex]=similarity
            table.loc[database[index][2].GetProp('_Name'),database[jndex][2].GetProp('_Name')]=similarity
    return table,hmap

smiles_char = ['?', '#', '%', ')', '(', '[', ']', '+', '-', '.', '0', '1', '2', '3', '4', '5',
               '6', '7', '8', '9', '=', 'A', 'C', 'B', 'E', 'D', 'G', 'F', 'I',
               'H', 'K', 'M', 'L', 'O', 'N', 'P', 'S', 'R', 'U', 'T', 'W', 'V',
               'Y', '[', 'Z', ']', '_', 'a', 'c', 'b', 'e', 'd', 'g', 'f', 'i',
               'h', 'm', 'l', 'o', 'n', 's', 'r', 'u', 't', 'y']  # 65

nucleic_char = ['A', 'U', 'C', 'G', 'N'] #5
structure_char = ['.', '(', ')', '[', ']', '?'] #6

from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
#enc_drug = OneHotEncoder().fit(np.array(smiles_char).reshape(-1, 1))
enc_drug = LabelEncoder().fit(np.array(smiles_char))
enc_protein = OneHotEncoder().fit(np.array(nucleic_char).reshape(-1, 1))
enc_structure = OneHotEncoder().fit(np.array(structure_char).reshape(-1, 1))

def drug_2_embed(x):
    #return enc_drug.transform(np.array(x).reshape(-1, 1)).toarray().T
    return enc_drug.transform(np.array(x))

def protein_2_embed(x):
    return enc_protein.transform(np.array(x).reshape(-1, 1)).toarray().T

def structure_2_embed(x):
    return enc_structure.transform(np.array(x).reshape(-1, 1)).toarray().T

def cluster_heatmap(database,table,hmap):
    '''
    进行层次聚类
    '''
    linked = linkage(hmap,'single')
    labelList = [mol[2].GetProp('_Name') for mol in database]

    plt.figure(figsize=(8,15),dpi=300)
    ax1=plt.subplot()
    o=dendrogram(linked,  
                orientation='left',
                labels=labelList,
                distance_sort='descending',
                show_leaf_counts=True)
    
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    plt.title('Similarity clustering',fontsize=20,weight='bold')
    plt.tick_params ('both',width=2,labelsize=8)
    plt.tight_layout()
    plt.savefig('docking/cluster_res/cluster_hierarchy.png') 


    # the clusters in order as the last plot
    new_data=list(reversed(o['ivl']))
    
    # create a new table with the order of HCL
    size=len(database)
    hmap_2=np.empty(shape=(size,size))
    for index,i in enumerate(new_data):
        for jndex,j in enumerate(new_data):
            hmap_2[index,jndex]=table.loc[i].at[j]
    
    figure= plt.figure(figsize=(30,30),dpi=300)
    gs1 = gridspec.GridSpec(2,7)
    gs1.update(wspace=0.01)
    ax1 = plt.subplot(gs1[0:-1, :2])
    dendrogram(linked, orientation='left', distance_sort='descending',show_leaf_counts=True,no_labels=True)
    ax1.spines['left'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    ax2 = plt.subplot(gs1[0:-1,2:6])
    f=ax2.imshow (hmap_2, cmap='PRGn_r', interpolation='nearest')
    
    ax2.set_title('Fingerprint Similarity',fontsize=20,weight='bold')
    ax2.set_xticks (range(len(new_data)))
    ax2.set_yticks (range(len(new_data)))
    ax2.set_xticklabels (new_data,rotation=90,size=8)
    ax2.set_yticklabels (new_data,size=8)
    
    ax3 = plt.subplot(gs1[0:-1,6:7])
    m=plt.colorbar(f,cax=ax3,shrink=0.75,orientation='vertical',spacing='uniform',pad=0.01)
    m.set_label ('Fingerprint Similarity')
    
    plt.tick_params ('both',width=2)
    plt.savefig('docking/cluster_res/cluster_heatmap.png')

def cluster_PCA_tSNE(database,flag='drug'):
    '''
    采用降维进行聚类
    '''
    target_value = list(set([1,0]))
    color_dict = {}
    colors = ['blue', 'red', 'yellow', 'green', 'orange', 'black', 'magenta', 'slategray', 'cyan', 'aquamarine']
    for i, t in enumerate(target_value):
        color_dict[t] = colors[i]
    print(color_dict)

    if flag=='drug':
        embeddings = []
        target = []
        for mol in database:
            #embeddings.append(drug_2_embed(list(mol[1])))
            embeddings.append(np.array(mol[3]))
            target.append(color_dict[mol[4]])
        embeddings = np.array(embeddings)
        target = np.array(target)
    elif flag=='RNA':
        embeddings = []
        target = []
        for mol in database:
            embeddings.append(np.array(protein_2_embed(list(mol[1]))))
            target.append(color_dict[int(mol[-1])])
        embeddings = np.array(embeddings)
        embeddings = np.reshape(embeddings,[embeddings.shape[0],embeddings.shape[1]*embeddings.shape[2]])
        target = np.array(target)
    else:
        return
    
    X_tsne = TSNE(learning_rate=200.0, perplexity=50).fit_transform(embeddings)
    X_pca = PCA().fit_transform(embeddings)

    plt.figure(figsize=(10, 5),dpi=300)
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1],s=3,c=target)
    plt.title('t-SNE')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('docking/cluster_res/cluster_tSNE.png')
    plt.close()

    plt.figure(figsize=(10, 5),dpi=300)
    plt.subplot(121)
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1],s=3,c=target)
    plt.title('t-SNE')
    plt.xticks([])
    plt.yticks([])

    plt.subplot(122)
    plt.scatter(X_pca[:, 0], X_pca[:, 1],s=3,c=target)
    plt.title('PCA')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('docking/cluster_res/cluster_PCA.png')

if __name__ == '__main__':
    print('run!')
    #prepareData('RNA')
    #calcDistribution('ALL')
    #database = loadingDrug('docking/data/ALL_valid_1.txt')
    database = loadingRNA('docking/data/ALL_valid_3.txt')
    #table,hmap = calcMatrix(database)
    #cluster_heatmap(database,table,hmap)
    cluster_PCA_tSNE(database,'RNA')
    print('end!')



