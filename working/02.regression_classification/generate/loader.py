"""
# Author: Yuhan Fei
# Created Time :  15 May 2020
# Revised Time :  25 May 2020
# File Name: loader.py
# Description: Data loader
"""

import pandas as pd
import numpy as np
import random
from torch.utils import data
from rdkit import Chem



def valid(input='./data/test.txt', output='./data/revised.txt'):
    print('This function is used to remove invalid SMILES...')
    f = open(output, "w")
    total=0
    removed=0
    for line in open(input):
        total += 1
        elements = line.strip().split("\t")
        if Chem.MolFromSmiles(elements[0]) != None:
            #drug = Chem.MolToSmiles(Chem.MolFromSmiles(elements[0]), isomericSmiles=isomeric)
            #f.write(str(drug) + "\t" + str(elements[1]) + "\n")
            f.write(str(elements[0]) + "\t" + str(elements[1]) + "\n")
        else:
            removed += 1
    print('Load %d records, remove %d invalid smiles...' % (total, removed))
    print('Remain %d records finally...' % (total-removed))
    print('Completed!')
    f.close()





def valid_ss(input='./data/test.txt', output='./data/revised.txt'):
    print('This function is used to remove invalid SMILES...')
    f = open(output, "w")
    total=0
    removed=0
    for line in open(input):
        total += 1
        elements = line.strip().split("\t")
        if Chem.MolFromSmiles(elements[0]) != None:
            #drug = Chem.MolToSmiles(Chem.MolFromSmiles(elements[0]), isomericSmiles=isomeric)
            #f.write(str(drug) + "\t" + str(elements[1]) + "\n")
            f.write(str(elements[0]) + "\t" + str(elements[1]) + "\t" + str(elements[2])+ "\t" + str(elements[3]) + "\n")
        else:
            removed += 1
    print('Load %d records, remove %d invalid smiles...' % (total, removed))
    print('Remain %d records finally...' % (total-removed))
    print('Completed!')
    f.close()


def valid_ss_knot(input='./data/test.txt', output='./data/revised.txt'):
    print('This function is used to remove invalid SMILES...')
    f = open(output, "w")
    total=0
    removed=0
    for line in open(input):
        total += 1
        elements = line.strip().split("\t")
        if Chem.MolFromSmiles(elements[0]) != None:
            #drug = Chem.MolToSmiles(Chem.MolFromSmiles(elements[0]), isomericSmiles=isomeric)
            #f.write(str(drug) + "\t" + str(elements[1]) + "\n")
            f.write(str(elements[0]) + "\t" + str(elements[1]) + "@" + str(elements[2]) +"\n")
        else:
            removed += 1
    print('Load %d records, remove %d invalid smiles...' % (total, removed))
    print('Remain %d records finally...' % (total-removed))
    print('Completed!')
    f.close()





def df2tab(input='./data/interaction.txt', output='./data/training.txt'):
    print('This function is used to convert two columns tabular format to complete adjacency matrix...')
    f = open(output, "w")
    df = pd.read_csv(input, sep='\t', header=None)
    total = len(df)
    first_col = len(df.take([0], axis=1).drop_duplicates())
    second_col = len(df.take([1], axis=1).drop_duplicates())
    df.columns = ['SMILES', 'SEQ']
    df2tab = lambda data: pd.get_dummies(data.SEQ).T.dot(data.SMILES.str.get_dummies(','))
    tab = df2tab(df)
    print('Load %d records,which contain %d SMILES and %d RNAs...' % (total, first_col, second_col))
    print('Convert complete adjacency matrix to complete two columns tabular format...')
    for index1, row in tab.iterrows():
        for index2, col in tab.iteritems():
            f.write(str(index2) + "\t" + str(index1) + "\t" + str(tab.at[index1, index2]) + "\n")
    f.close()
    df2 = pd.read_csv(output, sep='\t', header=None)
    total2 = len(df2)
    first_col2 = len(df2.take([0], axis=1).drop_duplicates())
    second_col2 = len(df2.take([1], axis=1).drop_duplicates())
    print('Generate %d records,which contain %d SMILES and %d RNAs ...' % (total2, first_col2, second_col2))
    print('Completed!')








def train_run(tab,output):
    f = open(output, "w")
    for index1, row in tab.iterrows():
        for index2, col in tab.iteritems():
            f.write(str(index1) + "\t" + str(index2) + "\t" + str(tab.at[index1, index2]) + "\n")
    f.close()




def df2var_rad(tab,weight=1):
    print('This function is used to generate dataset from tabular file...')
    print('The total number between positive and negative data is 1:%d' % (weight))
    X_drug_pos = []
    X_target_pos = []
    y_pos = []
    X_drug_neg = []
    X_target_neg = []
    y_neg = []

    for index1, row in tab.iterrows():
        for index2, col in tab.iteritems():
            if int(tab.at[index1, index2]) == 1:
                X_drug_pos.append(index1)
                X_target_pos.append(index2)
                y_pos.append(int(tab.at[index1, index2]))
            elif int(tab.at[index1, index2]) == 0:
                X_drug_neg.append(index1)
                X_target_neg.append(index2)
                y_neg.append(int(tab.at[index1, index2]))

    neg_tuple = list(zip(X_drug_neg, X_target_neg, y_neg))
    neg_tuple_random = random.sample(neg_tuple, int(len(X_drug_pos) * weight))
    X_drug_neg_final, X_target_neg_final, y_neg_final = map(list, zip(*neg_tuple_random))
    X_drug_final = X_drug_pos + X_drug_neg_final
    X_target_final = X_target_pos + X_target_neg_final
    y_final = y_pos + y_neg_final
    X_drug_num = len(X_drug_pos)
    X_drug_neg_num = len(X_drug_neg_final)
    y_num = len(y_final)
    print('Generate %d records,which contain %d postive records and %d negative records...' % (
    y_num, X_drug_num, X_drug_neg_num))
    print('Completed!')
    return np.array(X_drug_final), np.array(X_target_final), np.array(y_final)


def df2var_bal(tab,weight=1, random_seed=1):
    print('This function is used to generate dataset from tabular file...')
    print('The total number between positive and negative data is 1:%d' % (weight))
    X_drug_pos = []
    X_target_pos = []
    y_pos = []
    X_drug_neg = []
    X_target_neg = []
    y_neg = []

    for index1, row in tab.iterrows():
        for index2, col in tab.iteritems():
            if int(tab.at[index1, index2]) == 1:
                X_drug_pos.append(index1)
                X_target_pos.append(index2)
                y_pos.append(int(tab.at[index1, index2]))
            elif int(tab.at[index1, index2]) == 0:
                X_drug_neg.append(index1)
                X_target_neg.append(index2)
                y_neg.append(int(tab.at[index1, index2]))

    neg_tuple = list(zip(X_drug_neg, X_target_neg, y_neg))

    freq = {}
    for item in X_drug_pos:
        if (item in freq):
            freq[item] += 1
        else:
            freq[item] = 1
    X_drug_pos_freq = []
    for key, value in freq.items():
        a = [key, value]
        X_drug_pos_freq.append(a)


    neg_tuple_random = []
    for i in range(len(X_drug_pos_freq)):
        neg_tmp = []
        random.seed(random_seed+i)
        for j in range(len(neg_tuple)):
            if (X_drug_pos_freq[i][0] == neg_tuple[j][0]):
                neg_tmp.append(neg_tuple[j])

        neg_tuple_random.extend(random.sample(neg_tmp, X_drug_pos_freq[i][1] * weight))
    X_drug_neg_final, X_target_neg_final, y_neg_final = map(list, zip(*neg_tuple_random))




    X_drug_final = X_drug_pos + X_drug_neg_final
    X_target_final = X_target_pos + X_target_neg_final
    y_final = y_pos + y_neg_final
    X_drug_num = len(X_drug_pos)
    X_drug_neg_num = len(X_drug_neg_final)
    y_num = len(y_final)
    print('Generate %d records,which contain %d postive records and %d negative records...' % (
    y_num, X_drug_num, X_drug_neg_num))
    print('Completed!')
    return np.array(X_drug_final), np.array(X_target_final), np.array(y_final)







def df2var_bal_ss(tab,weight=1, random_seed=1):
    print('This function is used to generate dataset from tabular file...')
    print('The total number between positive and negative data is 1:%d' % (weight))
    X_drug_pos = []
    X_target_pos = []
    X_structure_pos =[]
    y_pos = []
    X_drug_neg = []
    X_target_neg = []
    X_structure_neg = []
    y_neg = []

    for index1, row in tab.iterrows():
        for index2, col in tab.iteritems():
            if int(tab.at[index1, index2]) == 1:
                X_drug_pos.append(index1)
                X_target_pos.append(index2.split("@")[0])
                X_structure_pos.append(index2.split("@")[1])
                y_pos.append(int(tab.at[index1, index2]))
            elif int(tab.at[index1, index2]) == 0:
                X_drug_neg.append(index1)
                X_target_neg.append(index2.split("@")[0])
                X_structure_neg.append(index2.split("@")[1])
                y_neg.append(int(tab.at[index1, index2]))

    neg_tuple = list(zip(X_drug_neg, X_target_neg, X_structure_neg, y_neg))

    freq = {}
    for item in X_drug_pos:
        if (item in freq):
            freq[item] += 1
        else:
            freq[item] = 1
    X_drug_pos_freq = []
    for key, value in freq.items():
        a = [key, value]
        X_drug_pos_freq.append(a)


    neg_tuple_random = []
    for i in range(len(X_drug_pos_freq)):
        neg_tmp = []
        random.seed(random_seed+i)
        #random.seed(random_seed)
        for j in range(len(neg_tuple)):
            if (X_drug_pos_freq[i][0] == neg_tuple[j][0]):
                neg_tmp.append(neg_tuple[j])
        neg_tuple_random.extend(random.sample(neg_tmp, X_drug_pos_freq[i][1] * weight))
    X_drug_neg_final, X_target_neg_final, X_structure_neg_final, y_neg_final = map(list, zip(*neg_tuple_random))




    X_drug_final = X_drug_pos + X_drug_neg_final
    X_target_final = X_target_pos + X_target_neg_final
    X_structure_final = X_structure_pos + X_structure_neg_final
    y_final = y_pos + y_neg_final
    X_drug_num = len(X_drug_pos)
    X_drug_neg_num = len(X_drug_neg_final)
    y_num = len(y_final)
    print('Generate %d records,which contain %d postive records and %d negative records...' % (
    y_num, X_drug_num, X_drug_neg_num))
    print('Completed!')
    return np.array(X_drug_final), np.array(X_target_final), np.array(X_structure_final), np.array(y_final)






def df2var_bal_ss2(tab,weight=1, random_seed=1):
    print('This function is used to generate dataset from tabular file...')
    print('The total number between positive and negative data is 1:%d' % (weight))
    X_drug_pos = []
    X_target_pos = []
    X_structure_pos =[]
    X_knot_pos =[]
    y_pos = []
    X_drug_neg = []
    X_target_neg = []
    X_structure_neg = []
    X_knot_neg = []
    y_neg = []

    for index1, row in tab.iterrows():
        for index2, col in tab.iteritems():
            if int(tab.at[index1, index2]) == 1:
                X_drug_pos.append(index1)
                X_target_pos.append(index2.split("@")[0])
                X_structure_pos.append(index2.split("@")[1])
                X_knot_pos.append(index2.split("@")[2])
                y_pos.append(int(tab.at[index1, index2]))
            elif int(tab.at[index1, index2]) == 0:
                X_drug_neg.append(index1)
                X_target_neg.append(index2.split("@")[0])
                X_structure_neg.append(index2.split("@")[1])
                X_knot_neg.append(index2.split("@")[2])
                y_neg.append(int(tab.at[index1, index2]))

    neg_tuple = list(zip(X_drug_neg, X_target_neg, X_structure_neg, X_knot_neg, y_neg))

    freq = {}
    for item in X_drug_pos:
        if (item in freq):
            freq[item] += 1
        else:
            freq[item] = 1
    X_drug_pos_freq = []
    for key, value in freq.items():
        a = [key, value]
        X_drug_pos_freq.append(a)


    neg_tuple_random = []
    for i in range(len(X_drug_pos_freq)):
        neg_tmp = []
        random.seed(random_seed+i)
        #random.seed(random_seed)
        for j in range(len(neg_tuple)):
            if (X_drug_pos_freq[i][0] == neg_tuple[j][0]):
                neg_tmp.append(neg_tuple[j])
        neg_tuple_random.extend(random.sample(neg_tmp, X_drug_pos_freq[i][1] * weight))
    X_drug_neg_final, X_target_neg_final, X_structure_neg_final, X_knot_neg_final, y_neg_final = map(list, zip(*neg_tuple_random))




    X_drug_final = X_drug_pos + X_drug_neg_final
    X_target_final = X_target_pos + X_target_neg_final
    X_structure_final = X_structure_pos + X_structure_neg_final
    X_knot_final = X_knot_pos + X_knot_neg_final
    y_final = y_pos + y_neg_final
    X_drug_num = len(X_drug_pos)
    X_drug_neg_num = len(X_drug_neg_final)
    y_num = len(y_final)
    print('Generate %d records,which contain %d postive records and %d negative records...' % (
    y_num, X_drug_num, X_drug_neg_num))
    print('Completed!')
    return np.array(X_drug_final), np.array(X_target_final), np.array(X_structure_final), np.array(X_knot_final), np.array(y_final)







def file2var_split(input, weight=1):
    print('This function is used to generate dataset from tabular file...')
    print('The total number between positive and negative data is 1:%d' % (weight))
    X_drug_pos = []
    X_target_pos = []
    y_pos = []
    X_drug_neg = []
    X_target_neg = []
    y_neg = []

    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        if int(values[2]) == 1:
            X_drug_pos.append(values[0])
            X_target_pos.append(values[1])
            y_pos.append(float(values[2]))
        elif int(values[2]) == 0:
            X_drug_neg.append(values[0])
            X_target_neg.append(values[1])
            y_neg.append(int(values[2]))
    file.close()
    neg_tuple = list(zip(X_drug_neg, X_target_neg, y_neg))
    neg_tuple_random = random.sample(neg_tuple, int(len(X_drug_pos) * weight))
    X_drug_neg_final, X_target_neg_final, y_neg_final = map(list, zip(*neg_tuple_random))
    X_drug_final = X_drug_pos + X_drug_neg_final
    X_target_final = X_target_pos + X_target_neg_final
    y_final = y_pos + y_neg_final
    X_drug_num = len(X_drug_pos)
    X_drug_neg_num = len(X_drug_neg_final)
    y_num = len(y_final)
    print('Generate %d records,which contain %d postive records and %d negative records...' % (
    y_num, X_drug_num, X_drug_neg_num))
    print('Completed!')
    return np.array(X_drug_final), np.array(X_target_final), np.array(y_final)


def file2var_rad(input, weight=1):
    print('This function is used to generate dataset from tabular file...')
    print('The total number between positive and negative data is 1:%d' % (weight))
    X_drug_pos = []
    X_target_pos = []
    y_pos = []
    X_drug_neg = []
    X_target_neg = []
    y_neg = []

    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        if int(values[2]) == 1:
            X_drug_pos.append(values[0])
            X_target_pos.append(values[1])
            y_pos.append(float(values[2]))
        elif int(values[2]) == 0:
            X_drug_neg.append(values[0])
            X_target_neg.append(values[1])
            y_neg.append(int(values[2]))
    file.close()
    neg_tuple = list(zip(X_drug_neg, X_target_neg, y_neg))
    neg_tuple_random = random.sample(neg_tuple, int(len(X_drug_pos) * weight))
    X_drug_neg_final, X_target_neg_final, y_neg_final = map(list, zip(*neg_tuple_random))
    X_drug_final = X_drug_pos + X_drug_neg_final
    X_target_final = X_target_pos + X_target_neg_final
    y_final = y_pos + y_neg_final
    X_drug_num = len(X_drug_pos)
    X_drug_neg_num = len(X_drug_neg_final)
    y_num = len(y_final)
    print('Generate %d records,which contain %d postive records and %d negative records...' % (
    y_num, X_drug_num, X_drug_neg_num))
    print('Completed!')
    return np.array(X_drug_final), np.array(X_target_final), np.array(y_final)




def file2var_bal(input, weight=1):
    print('This function is used to generate dataset from tabular file with same degree...')
    print('The total number between positive and negative data is 1:%d' % (weight))
    df = pd.read_table(input, delimiter = "\t",header=None)
    df.columns = ['SMILES', 'Seq','Label']
    df1=df.loc[df['Label'] == 1]
    df2=df1['SMILES'].value_counts()
    df3=df2.to_frame().reset_index()
    df3.columns = ['SMILES', 'Freq']
    column_names = ["SMILES", "Seq", "Label"]
    neg = pd.DataFrame(columns = column_names, dtype=object)
    for row in df3.itertuples():
        SMILES=getattr(row, 'SMILES')
        Freq=getattr(row, 'Freq')*weight
        a=df.loc[(df['SMILES'] == SMILES) & (df['Label'] ==0)]
        b=a.sample(n=Freq, frac=None, replace=False, weights=None, random_state=None, axis=0)
        neg=neg.append(b,ignore_index = True)
    frames = [df1, neg]
    result = pd.concat(frames)
    print('Generate %d records,which contain %d postive records and %d negative records...' % (
    result.shape[0], df1.shape[0], neg.shape[0]))
    print('Completed!')
    return result['SMILES'].values, result['Seq'].values, result['Label'].values



def file2var(input):
    print('This function is used to generate dataset from tabular file...')
    X_drug = []
    X_target = []
    y = []

    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        X_drug.append(values[0])
        X_target.append(values[1])
        y.append(float(values[2]))
    file.close()
    print('Completed!')
    return np.array(X_drug), np.array(X_target), np.array(y)


def convert_y_unit(y, from_, to_):
    array_flag = False
    if isinstance(y, (int, float)):
        y = np.array([y])
        array_flag = True
    y = y.astype(float)    
    # basis as nM
    if from_ == 'nM':
        y = y
    elif from_ == 'p':
        y = 10**(-y) / 1e-9

    if to_ == 'p':
        zero_idxs = np.where(y == 0.)[0]
        y[zero_idxs] = 1e-10
        y = -np.log10(y*1e-9)
    elif to_ == 'nM':
        y = y
        
    if array_flag:
        return y[0]
    return y


def file2var_reg(input, convert_to_log = True, ):
    print('This function is used to generate dataset from tabular file...')
    X_drug = []
    X_target = []
    y = []
    if convert_to_log:
        print('Default set to logspace (nM -> p) for easier regression')


    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        X_drug.append(values[0])
        X_target.append(values[1])
        if convert_to_log:
            y.append(convert_y_unit(np.array(float(values[2])), 'nM', 'p'))
        else:
            y.append(float(values[2]))
    file.close()
    print('Completed!')
    return np.array(X_drug), np.array(X_target), np.array(y)





def file2var_reg_ss(input, convert_to_log = True, ):
    print('This function is used to generate dataset from tabular file...')
    X_drug = []
    X_target = []
    X_structure = []    
    y = []
    if convert_to_log:
        print('Default set to logspace (nM -> p) for easier regression')


    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        X_drug.append(values[0])
        X_target.append(values[1])
        X_structure.append(values[2])
        if convert_to_log:
            y.append(convert_y_unit(np.array(float(values[3])), 'nM', 'p'))
        else:
            y.append(float(values[3]))
    file.close()
    print('Completed!')
    return np.array(X_drug), np.array(X_target), np.array(X_structure), np.array(y)


def file2var_reg_ss2(input, convert_to_log = True, ):
    print('This function is used to generate dataset from tabular file...')
    X_drug = []
    X_target = []
    X_structure = []
    X_knot = []     
    y = []
    if convert_to_log:
        print('Default set to logspace (nM -> p) for easier regression')


    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        X_drug.append(values[0])
        X_target.append(values[1])
        X_structure.append(values[2])
        X_knot.append(values[3])
        if convert_to_log:
            y.append(convert_y_unit(np.array(float(values[4])), 'nM', 'p'))
        else:
            y.append(float(values[4]))
    file.close()
    print('Completed!')
    return np.array(X_drug), np.array(X_target), np.array(X_structure), np.array(X_knot), np.array(y)




def file2var_reg_ss_aff(input, convert_to_log = True):
    print('This function is used to generate dataset from tabular file...')
    X_drug = []
    X_target = []
    X_structure = []  
    X_affinity = []  
    y = []
    if convert_to_log:
        print('Default set to logspace (nM -> p) for easier regression')


    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        X_drug.append(values[0])
        X_target.append(values[1])
        X_structure.append(values[2])
        X_affinity.append(values[3])
        if convert_to_log:
            y.append(convert_y_unit(np.array(float(values[4])), 'nM', 'p'))
        else:
            y.append(float(values[4]))
    file.close()
    print('Completed!')
    return np.array(X_drug), np.array(X_target), np.array(X_structure), np.array(X_affinity), np.array(y)





def file2var_ss(input):
    print('This function is used to generate dataset from tabular file...')
    X_drug = []
    X_target = []
    X_structure = []
    y = []

    file = open(input, "r")
    for aline in file:
        values = aline.split("\t")
        X_drug.append(values[0])
        X_target.append(values[1])
        X_structure.append(values[2])
        y.append(float(values[3]))
    file.close()
    print('Completed!')
    return np.array(X_drug), np.array(X_target), np.array(X_structure), np.array(y)



def split(input='./data/interaction.txt', frac=[0.7,0.1,0.2],random_state=1):
    print('This function is used to convert two columns tabular format to complete adjacency matrix...')
    df = pd.read_csv(input, sep='\t', header=None)
    df.columns = ['SMILES', 'SEQ']

    print('splitting dataset...')
    train_frac, val_frac, test_frac = frac
    test = df.sample(frac=test_frac, replace=False, random_state=random_state)
    train_val = df[~df.index.isin(test.index)]
    val = train_val.sample(frac=val_frac / (1 - test_frac), replace=False, random_state=random_state)
    train = train_val[~train_val.index.isin(val.index)]


    f1 = open("./train_bal.txt", "w")
    df2tab = lambda data: pd.get_dummies(data.SEQ).T.dot(data.SMILES.str.get_dummies(','))
    train_tab = df2tab(train)
    for index1, row in train_tab.iterrows():
        for index2, col in train_tab.iteritems():
            f1.write(str(index2) + "\t" + str(index1) + "\t" + str(train_tab.at[index1, index2]) + "\n")
    f1.close()
    print('Completed_train!')
    f2 = open("./val_bal.txt", "w")
    df2tab = lambda data: pd.get_dummies(data.SEQ).T.dot(data.SMILES.str.get_dummies(','))
    val_tab = df2tab(val)
    for index1, row in val_tab.iterrows():
        for index2, col in val_tab.iteritems():
            f2.write(str(index2) + "\t" + str(index1) + "\t" + str(val_tab.at[index1, index2]) + "\n")
    f2.close()
    print('Completed_validation!')
    f3 = open("./test_bal.txt", "w")
    df2tab = lambda data: pd.get_dummies(data.SEQ).T.dot(data.SMILES.str.get_dummies(','))
    test_tab = df2tab(test)
    for index1, row in test_tab.iterrows():
        for index2, col in test_tab.iteritems():
            f3.write(str(index2) + "\t" + str(index1) + "\t" + str(test_tab.at[index1, index2]) + "\n")
    f3.close()
    print('Completed_test!')








def df2tab(input='./data/interaction.txt', output='./data/training.txt'):
    print('This function is used to convert two columns tabular format to complete adjacency matrix...')
    f = open(output, "w")
    df = pd.read_csv(input, sep='\t', header=None)
    total = len(df)
    first_col = len(df.take([0], axis=1).drop_duplicates())
    second_col = len(df.take([1], axis=1).drop_duplicates())
    df.columns = ['SMILES', 'SEQ']
    df2tab = lambda data: pd.get_dummies(data.SEQ).T.dot(data.SMILES.str.get_dummies(','))
    tab = df2tab(df)
    print('Load %d records,which contain %d SMILES and %d RNAs...' % (total, first_col, second_col))
    print('Convert complete adjacency matrix to complete two columns tabular format...')
    for index1, row in tab.iterrows():
        for index2, col in tab.iteritems():
            f.write(str(index2) + "\t" + str(index1) + "\t" + str(tab.at[index1, index2]) + "\n")
    f.close()
    df2 = pd.read_csv(output, sep='\t', header=None)
    total2 = len(df2)
    first_col2 = len(df2.take([0], axis=1).drop_duplicates())
    second_col2 = len(df2.take([1], axis=1).drop_duplicates())
    print('Generate %d records,which contain %d SMILES and %d RNAs ...' % (total2, first_col2, second_col2))
    print('Completed!')










def split_drug(input='./data/interaction.txt', frac=[0.7,0.1,0.2],random_seed=1):
    print('This function is used to convert two columns tabular format to complete adjacency matrix...')
    df = pd.read_csv(input, sep='\t', header=None)
    df.columns = ['SMILES', 'SEQ']
    df2tab = lambda data: pd.get_dummies(data.SMILES).T.dot(data.SEQ.str.get_dummies(','))
    tab = df2tab(df)



    print('splitting dataset...')
    train_frac, val_frac, test_frac = frac
    test = tab.sample(frac=test_frac, replace=False, random_state=random_seed)
    train_val = tab[~tab.index.isin(test.index)]
    val = train_val.sample(frac=val_frac / (1 - test_frac), replace=False, random_state=random_seed)
    train = train_val[~train_val.index.isin(val.index)]

    return train, val, test




def split_target(input='./data/interaction.txt', frac=[0.7,0.1,0.2],random_state=1):
    print('This function is used to convert two columns tabular format to complete adjacency matrix...')
    df = pd.read_csv(input, sep='\t', header=None)
    df.columns = ['SMILES', 'SEQ']
    df2tab = lambda data: pd.get_dummies(data.SEQ).T.dot(data.SMILES.str.get_dummies(','))
    tab = df2tab(df)


    print('splitting dataset...')
    train_frac, val_frac, test_frac = frac
    test = tab.sample(frac=test_frac, replace=False, random_state=random_state)
    train_val = tab[~tab.index.isin(test.index)]
    val = train_val.sample(frac=val_frac / (1 - test_frac), replace=False, random_state=random_state)
    train = train_val[~train_val.index.isin(val.index)]


    return train, val, test


    """
    X_drugs_train, X_targets_train, y_train = loader.file2var_rad(input='train_run.txt', weight=2)
    X_drugs_val, X_targets_val, y_val = loader.file2var_rad(input='val_run.txt', weight=2)
    X_drugs_test, X_targets_test, y_test = loader.file2var_rad(input='test_run.txt', weight=2)
    """


#amino_char = ['?', 'A', 'C', 'B', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'O','N', 'Q', 'P', 'S', 'R', 'U', 'T', 'W', 'V', 'Y', 'X', 'Z']  # 26

#nucleic_char = ['A', 'U', 'C', 'G', 'X','N']

#nucleic_char = ['A', 'U', 'C', 'G', 'N']



smiles_char = ['?', '#', '%', ')', '(', '[', ']', '+', '-', '.', '0', '1', '2', '3', '4', '5',
               '6', '7', '8', '9', '=', 'A', 'C', 'B', 'E', 'D', 'G', 'F', 'I',
               'H', 'K', 'M', 'L', 'O', 'N', 'P', 'S', 'R', 'U', 'T', 'W', 'V',
               'Y', '[', 'Z', ']', '_', 'a', 'c', 'b', 'e', 'd', 'g', 'f', 'i',
               'h', 'm', 'l', 'o', 'n', 's', 'r', 'u', 't', 'y']  # 67

nucleic_char = ['A', 'U', 'C', 'G', 'N' ] #4


#structure_char = ['.', '(']  # 2

#structure_char = [ '.', '(', '[' ]  # 3

#structure_char = [ '.', '(', ')' ]  # 3

#structure_char = ['.', '(', '?']  # 3

structure_char = ['.', '(', '[', '?']  # 4

#structure_char = ['.', '(', 'J','L','H','I']  #6

#knot_char = [ '.', '[' ,'{','<','A','B','C','D','E']

knot_char = [ '.', '[']

def predict_load(path, n=50):
	# a line in the file is SMILES Target_seq score
	try:
		file = open(path, "r")
	except:
		print('Path Not Found, please double check!')
	X_drug = []
	X_target = []
	for line in file.readlines()[0:n]:
		values = line.split('\t')
		X_drug.append(values[0])
		X_target.append(values[1])
	file.close()
	return np.array(X_drug), np.array(X_target)


def predict_load_ss(path, n=50):
    # a line in the file is SMILES Target_seq score
    try:
        file = open(path, "r")
    except:
        print('Path Not Found, please double check!')
    X_drug = []
    X_target = []
    X_ss = []
    for line in file.readlines()[0:n]:
        values = line.split('\t')
        X_drug.append(values[0])
        X_target.append(values[1])
        X_ss.append(values[2])
    file.close()
    return np.array(X_drug), np.array(X_target), np.array(X_ss)


def predict_load_ss_aff(path):
    # a line in the file is SMILES Target_seq score
    try:
        file = open(path, "r")
    except:
        print('Path Not Found, please double check!')
    X_tname = []
    X_target = []
    X_ss = []
    X_dname = []
    X_drug = []
    X_aff = []
    for line in file.readlines():
        values = line.split('\t')
        X_tname.append(values[0])
        X_target.append(values[1])
        X_ss.append(values[2])
        X_dname.append(values[3])
        X_drug.append(values[4])
        X_aff.append(values[5].replace("\n",""))
    file.close()
    return np.array(X_tname), np.array(X_target), np.array(X_ss), np.array(X_dname), np.array(X_drug), np.array(X_aff)




from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder


enc_protein = OneHotEncoder().fit(np.array(nucleic_char).reshape(-1, 1))
enc_drug = OneHotEncoder().fit(np.array(smiles_char).reshape(-1, 1))
enc_structure = OneHotEncoder().fit(np.array(structure_char).reshape(-1, 1))

enc_structure2 =  LabelEncoder().fit(np.array(structure_char))




enc_knot2 =  LabelEncoder().fit(np.array(knot_char))


class data_process_loader(data.Dataset):

    def __init__(self, list_IDs, labels, df, **config):
        'Initialization'
        self.labels = labels
        self.list_IDs = list_IDs
        self.df = df
        self.config = config

        if self.config['drug_encoding'] in ['DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP']:
            from dgllife.utils import smiles_to_bigraph, CanonicalAtomFeaturizer, CanonicalBondFeaturizer
            self.node_featurizer = CanonicalAtomFeaturizer()
            self.edge_featurizer = CanonicalBondFeaturizer(self_loop=True)
            from functools import partial
            self.fc = partial(smiles_to_bigraph, add_self_loop=True)

        elif self.config['drug_encoding'] == 'DGL_AttentiveFP':
            from dgllife.utils import smiles_to_bigraph, AttentiveFPAtomFeaturizer, AttentiveFPBondFeaturizer
            self.node_featurizer = AttentiveFPAtomFeaturizer()
            self.edge_featurizer = AttentiveFPBondFeaturizer(self_loop=True)
            from functools import partial
            self.fc = partial(smiles_to_bigraph, add_self_loop=True)


        elif self.config['drug_encoding'] in ['DGL_GIN_AttrMasking', 'DGL_GIN_ContextPred']:
            from dgllife.utils import smiles_to_bigraph, PretrainAtomFeaturizer, PretrainBondFeaturizer
            self.node_featurizer = PretrainAtomFeaturizer()
            self.edge_featurizer = PretrainBondFeaturizer()
            from functools import partial
            self.fc = partial(smiles_to_bigraph, add_self_loop=True)

    def __len__(self):
        'Denotes the total number of samples'
        return len(self.list_IDs)

    def __getitem__(self, index):
        'Generates one sample of data'
        index = self.list_IDs[index]
        v_d = self.df.iloc[index]['drug_encoding']

        if self.config['drug_encoding'] == 'CNN' or self.config['drug_encoding'] == 'RNN':
            v_d = drug_2_embed(v_d)
        elif self.config['drug_encoding'] in ['DGL_GCN', 'DGL_GAT', 'DGL_NeuralFP', 'DGL_GIN_AttrMasking', 'DGL_GIN_ContextPred','DGL_AttentiveFP']:
            v_d = self.fc(smiles=v_d, node_featurizer=self.node_featurizer, edge_featurizer=self.edge_featurizer)
        v_p = self.df.iloc[index]['target_encoding']
        v_s = self.df.iloc[index]['structure_encoding']
        v_k = self.df.iloc[index]['knot_encoding']
        #v_a = self.df.iloc[index]['Affinity']
        #v_a=float(v_a)

        #v_a= convert_y_unit(np.array(float(v_a)), 'nM', 'p')

        #v_a=abs(float(v_a))

        #v_a = 1.0 / ( 1 + np.exp(-float(v_a)) )

        #v_a=np.asfarray([v_a]*31,float)
        #print(v_a)

        if self.config['target_encoding'] == 'CNN' or self.config['target_encoding'] == 'RNN' or self.config['target_encoding'] == 'Prism':
            #print(v_p)
            #print(v_s)
            #print(v_p)
            #print(v_p)
            v_p = protein_2_embed(v_p)
            v_s = structure_2_embed2(v_s)
            v_k = knot_2_embed2(v_k)

            #v_a = affinity_2_embed2(v_a)
        #print(v_s)
        #print(type(v_s))

        #v_p = np.expand_dims(v_p, axis=2).transpose([2, 1, 0])
        #v_s = np.expand_dims(v_s, axis=2).transpose([2, 1, 0])
        #v_ps = np.vstack((v_p,v_s))
        #v_ps = np.expand_dims(v_ps, axis=2).transpose([2, 1, 0])
        v_psk = np.vstack((v_p,v_s,v_k))
        #v_psk = np.expand_dims(v_psk, axis=2).transpose([2, 1, 0])
        #v_psa = np.vstack((v_p,v_s,v_a))

        #print(v_ps.shape)
        
       
        #print(v_pss.shape)
        y = self.labels[index]
        return v_d, v_psk, y
        #return v_d, v_p, y
        #return v_d, v_psa, y
        #return v_d, v_ps, y


def drug_2_embed(x):
    return enc_drug.transform(np.array(x).reshape(-1, 1)).toarray().T


def protein_2_embed(x):
    return enc_protein.transform(np.array(x).reshape(-1, 1)).toarray().T


def structure_2_embed(x):
    return enc_structure.transform(np.array(x).reshape(-1, 1)).toarray().T


def structure_2_embed2(x):
    return enc_structure2.transform(np.array(x))


def affinity_2_embed(x):
    return enc_affinity.transform(np.array(x).reshape(-1, 1)).toarray().T


def knot_2_embed2(x):
    return enc_knot2.transform(np.array(x))











def load_process(path = './data', binary = False, convert_to_log = True, threshold = 30):
    print('Beginning Processing...')



    print('Beginning to extract zip file...')
    with ZipFile(saved_path, 'r') as zip:
        zip.extractall(path = path)

    affinity = pd.read_csv(path + '/DAVIS/affinity.txt', header=None, sep = ' ')

    with open(path + '/DAVIS/target_seq.txt') as f:
        target = json.load(f)

    with open(path + '/DAVIS/SMILES.txt') as f:
        drug = json.load(f)

    target = list(target.values())
    drug = list(drug.values())

    SMILES = []
    Target_seq = []
    y = []

    for i in range(len(drug)):
        for j in range(len(target)):
            SMILES.append(drug[i])
            Target_seq.append(target[j])
            y.append(affinity.values[i, j])

    if binary:
        print('Default binary threshold for the binding affinity scores are 30, you can adjust it by using the "threshold" parameter')
        y = [1 if i else 0 for i in np.array(y) < threshold]
    else:
        if convert_to_log:
            print('Default set to logspace (nM -> p) for easier regression')
            y = convert_y_unit(np.array(y), 'nM', 'p')
        else:
            y = y
    print('Done!')
    return np.array(SMILES), np.array(Target_seq), np.array(y)
