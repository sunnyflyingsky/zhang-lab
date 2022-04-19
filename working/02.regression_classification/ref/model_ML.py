#!usr/bin/python
# -*- coding: utf-8 -*-
import sys, re
import numpy as np


if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parse = argparse.ArgumentParser(description='This is kmer module for generate nucleic acid compositio vector.', formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfiles', nargs='*', help='The input files in FASTA format.')
    parse.add_argument('alphabet', choices=['DNA', 'RNA'], help='The sequence type.')
    parse.add_argument('method', type=str, help='The method name of nucleic acid composition. {Kmer,mismatch,subsequence}')
    parse.add_argument('-k', type=int, help='For Kmer, IDKmer, mismatch, subsequence methods. The k value of kmer.')
    parse.add_argument('-r', type=int, choices=[1, 0], help='Whether consider the reverse complement or not.\n1 means True, 0 means False. (default = 0)')
    args = parse.parse_args()


    #print(main(args))

    Seqs, cats = main(args)
    #print(Seqs)


    X_resampled_smote = np.array(Seqs)




    from sklearn.preprocessing import LabelEncoder
    import warnings
    warnings.filterwarnings(module='sklearn*', action='ignore', category=DeprecationWarning)
    le = LabelEncoder()
    y_resampled_smote = le.fit_transform(cats)

    from collections import Counter
    print(sorted(Counter(y_resampled_smote).items()))


    AGO= np.count_nonzero(y_resampled_smote == 1)
    WT = np.count_nonzero(y_resampled_smote == 0)
    #All=int(AGO)+int(WT)
    if int(WT-AGO)>=0:
        Weight=float(WT/AGO)
    else:
        Weight=float(AGO/WT)


    from sklearn import model_selection
    Train_X, Test_X, Train_y, Test_y = model_selection.train_test_split(X_resampled_smote, y_resampled_smote, stratify=y_resampled_smote, test_size=0.2,random_state=2020)

    import numpy as np
    from sklearn.model_selection import cross_val_score


    from xgboost import XGBClassifier
    #from sklearn.svm import SVC
    #from sklearn.linear_model import LogisticRegression, SGDClassifier
    #from sklearn.svm import SVC
    #from sklearn.tree import DecisionTreeClassifier
    #from sklearn.ensemble import RandomForestClassifier
    #from sklearn.model_selection import cross_val_score
    #from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
    #from xgboost import XGBClassifier
    #import lightgbm as lgb #lgb.LGBMClassifier
    #from catboost import CatBoostClassifier
    #from sklearn.neighbors import KNeighborsClassifier


    def classifier():
        scores=[]
        f1s=[]
        aucs=[]
        precisions=[]
        recalls=[]
        model = XGBClassifier(random_state=2020, n_estimators=200, eta=0.1, max_depth=6, objective='binary:logistic',scale_pos_weight=float(Weight))
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


    scores, f1s, aucs, precisions, recalls = classifier()
    print('accurcay', scores)
    print('precision', precisions)
    print('recall', recalls)
    print('f1', f1s)




    #模型生成
    import pickle
    xgb = XGBClassifier(random_state=2020, n_estimators=200, eta=0.1, max_depth=6, objective='binary:logistic',scale_pos_weight=float(Weight))
    model=xgb.fit(Train_X, Train_y)
    f = open('model.pickle', 'wb')
    pickle.dump(model, f)
    f.close()

