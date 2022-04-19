#!usr/bin/python
# -*- coding: utf-8 -*-
import sys, re
import numpy as np
from sklearn import model_selection
from sklearn.preprocessing import LabelEncoder
from sklearn import preprocessing
import pandas as pd
from sklearn.model_selection import cross_val_score

from xgboost import XGBClassifier,XGBRegressor,XGBRFRegressor
import lightgbm as lgb #lgb.LGBMClassifier

#from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression,LinearRegression
from sklearn.ensemble import AdaBoostRegressor, GradientBoostingRegressor,BaggingRegressor,RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR

#from sklearn.svm import SVC
#from sklearn.tree import DecisionTreeClassifier
#from sklearn.model_selection import cross_val_score
#from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
#from catboost import CatBoostClassifier
#from sklearn.neighbors import KNeighborsClassifier

smiles_char = ['?', '#', '%', ')', '(', '[', ']', '+', '-', '.', '0', '1', '2', '3', '4', '5',
               '6', '7', '8', '9', '=', 'A', 'C', 'B', 'E', 'D', 'G', 'F', 'I',
               'H', 'K', 'M', 'L', 'O', 'N', 'P', 'S', 'R', 'U', 'T', 'W', 'V',
               'Y', '[', 'Z', ']', '_', 'a', 'c', 'b', 'e', 'd', 'g', 'f', 'i',
               'h', 'm', 'l', 'o', 'n', 's', 'r', 'u', 't', 'y']  # 65
nucleic_char = ['N','A', 'U', 'C', 'G'] #5
structure_char = ['?','.','[', '(']  # 4
#knot_char = [ '.', '[']

enc_smile = LabelEncoder().fit(np.array(smiles_char))
enc_nucleic = LabelEncoder().fit(np.array(nucleic_char))
enc_structure =  LabelEncoder().fit(np.array(structure_char))
#enc_knot =  LabelEncoder().fit(np.array(knot_char))

if __name__ == '__main__':
    #import argparse
    #from argparse import RawTextHelpFormatter
    #parse = argparse.ArgumentParser(description='This is kmer module for generate nucleic acid compositio vector.', formatter_class=RawTextHelpFormatter)
    #parse.add_argument('inputfiles', nargs='*', help='The input files in FASTA format.')
    #parse.add_argument('alphabet', choices=['DNA', 'RNA'], help='The sequence type.')
    #parse.add_argument('method', type=str, help='The method name of nucleic acid composition. {Kmer,mismatch,subsequence}')
    #arse.add_argument('-k', type=int, help='For Kmer, IDKmer, mismatch, subsequence methods. The k value of kmer.')
    #arse.add_argument('-r', type=int, choices=[1, 0], help='Whether consider the reverse complement or not.\n1 means True, 0 means False. (default = 0)')
    #args = parse.parse_args()
    #print(main(args))

    #Seqs, cats = main(args)
    #print(Seqs)
    #smartnet/smartnet_regression/data/final_data_autodock_knot.txt
    data = pd.read_csv('smartnet/smartnet_regression/data/final_data_dock6_knot.txt',sep='\t',names=['smiles','seqs','structs','labels'])
    Seqs,lables = [],[]
    for ix in data.index.values:
        #a = list(data.iloc[ix]['smiles'])
        s_d = data.iloc[ix]['smiles']
        s_r = data.iloc[ix]['seqs']
        s_k = data.iloc[ix]['structs']
        if len(s_d)<100:
            s_d += '?'*(100-len(s_d))
        else:
            s_d = s_d[:100]
        if len(s_r)<31:
            s_r += 'N'*(31-len(s_r))
        else:
            s_r = s_r[:31]
        if len(s_k)<31:
            s_k += 'N'*(31-len(s_k))
        else:
            s_k = s_k[:31]  
        v_d = enc_smile.transform(np.array(list(s_d))).tolist()
        v_r = enc_nucleic.transform(np.array(list(s_r))).tolist()
        v_k = enc_structure.transform(np.array(list(s_k))).tolist()
        #print(data.iloc[ix]['seqs'],v_r)
        v_l = data.iloc[ix]['labels']
        Seqs.append(v_d + v_r + v_k)
        #print(max(Seqs[-1]))
        lables.append(v_l)
    #Seqs = Seqs[:50]
    #lables = lables[:50]
    X_resampled_smote = np.array(Seqs)

    from sklearn.preprocessing import LabelEncoder
    import warnings
    warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
    #le = LabelEncoder()
    #y_resampled_smote = le.fit_transform(lables)
    y_resampled_smote = np.array(lables)

    """
    from collections import Counter
    print(sorted(Counter(y_resampled_smote).items()))


    AGO= np.count_nonzero(y_resampled_smote == 1)
    WT = np.count_nonzero(y_resampled_smote == 0)
    #All=int(AGO)+int(WT)
    if int(WT-AGO)>=0:
        Weight=float(WT/AGO)
    else:
        Weight=float(AGO/WT)
    """
    X_resampled_smote = preprocessing.scale(X_resampled_smote)
    #y_resampled_smote = preprocessing.scale(y_resampled_smote)
    Train_X, Test_X, Train_y, Test_y = model_selection.train_test_split(X_resampled_smote, y_resampled_smote, test_size=0.2,random_state=2020)

    def classifier():
        scores=[]
        f1s=[]
        aucs=[]
        precisions=[]
        recalls=[]
        model = XGBClassifier(random_state=2020, n_estimators=200, eta=0.1, max_depth=6, objective='reg:logistic')
        #model = SVC(class_weight='balanced',random_state=2020)

        score = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='accuracy')
        mean_score = np.round(np.mean(score),5)
        scores.append(mean_score)

        precision = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='precision')
        mean_precision = np.round(np.mean(precision),5)
        precisions.append(mean_precision)

        recall = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='recall')
        mean_recall = np.round(np.mean(recall),5)
        recalls.append(mean_recall)

        f1 = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='f1')
        mean_f1 = np.round(np.mean(f1),5)
        f1s.append(mean_f1)

        auc = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='roc_auc')
        mean_auc = np.round(np.mean(auc), 5)
        aucs.append(mean_auc)

        return scores, f1s, aucs, precisions, recalls
    
    #scores, f1s, aucs, precisions, recalls = classifier()
    #print('accurcay', scores)
    #print('precision', precisions)
    #print('recall', recalls)
    #print('f1', f1s)

    #模型生成
    import pickle
    xgb = XGBRFRegressor(random_state=2020, n_estimators=120, eta=0.1, max_depth=3, silent=False, objective='reg:gamma')  #multi:softprob, reg:logistic
    model=xgb.fit(Train_X, Train_y)
    #gbm = lgb.LGBMClassifier(num_leaves=31, learning_rate=0.05, n_estimators=20, objective='regression')
    #model = gbm.fit(Train_X, Train_y, eval_set=[(Test_X, Test_y)], eval_metric='l1', early_stopping_rounds=5)
    #LR=LogisticRegression()
    #model = LR.fit(Train_X,Train_y)
    #clf = RandomForestRegressor()
    #model = clf.fit(Train_X,Train_y)
    clf = GradientBoostingRegressor()
    model = clf.fit(Train_X,Train_y)
    print("Model Score_train -> ", model.score(Train_X, Train_y))
    print("Model Score_test -> ", model.score(Test_X, Test_y))
    f = open('smartnet/smartnet_vscode/result/xgb.pickle', 'wb')
    pickle.dump(model, f)
    f.close()


    from sklearn.model_selection import learning_curve
    import matplotlib.pyplot as plt
    #plt.switch_backend('agg')
    # 学习曲线

    def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None,
                            n_jobs=1, train_sizes=np.linspace(.1, 1.0, 5)):
        plt.figure()
        plt.title(title)
        if ylim is not None:
            plt.ylim(*ylim)
        plt.xlabel("Training examples")
        plt.ylabel("Score")
        train_sizes, train_scores, test_scores = learning_curve(estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
        train_scores_mean = np.mean(train_scores, axis=1)
        train_scores_std = np.std(train_scores, axis=1)
        test_scores_mean = np.mean(test_scores, axis=1)
        test_scores_std = np.std(test_scores, axis=1)
        plt.grid()

        plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                         train_scores_mean + train_scores_std, alpha=0.1,
                         color="r")
        plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                         test_scores_mean + test_scores_std, alpha=0.1, color="g")
        plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
                 label="Training score")
        plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
                 label="Cross-validation score")

        plt.legend(loc="lower left") #best
        return plt

    from sklearn.model_selection import ShuffleSplit

    title = "Learning Curves " #(SVM, RBF kernel, $\gamma=0.001$)
    cv = ShuffleSplit(n_splits=10, test_size=0.2, random_state=2020)
    # estimator = XGBClassifier(silent=True)  # 建模

    plot_learning_curve(xgb, title, Train_X, Train_y, (-2.0, 0.0), cv=cv, n_jobs=1)
    plt.savefig('smartnet/smartnet_vscode/result/lc_xgb.png')
