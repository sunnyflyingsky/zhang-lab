

if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parse = argparse.ArgumentParser(description='This is kmer module for generate nucleic acid compositio vector.', formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfiles', nargs='*', help='The input files in FASTA format.')
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'], help='The sequence type.')
    parse.add_argument('method', type=str, help='The method name of nucleic acid composition. {Kmer,mismatch,subsequence}')
    parse.add_argument('-k', type=int, help='For Kmer, IDKmer, mismatch, subsequence methods. The k value of kmer.')
    parse.add_argument('-r', type=int, choices=[1, 0], help='Whether consider the reverse complement or not.\n1 means True, 0 means False. (default = 0)')
    args = parse.parse_args()


    #print(main(args))

    Seqs, cats = main(args)
    #print(Seqs)



    X_resampled_smote = np.array(Seqs)


    #均一化
    #from sklearn.preprocessing import StandardScaler
    #scaler = StandardScaler()
    #X_resampled = scaler.fit_transform(X_resampled)

    #from sklearn.preprocessing import MinMaxScaler
    #scaler = MinMaxScaler()
    #X_resampled = scaler.fit_transform(X_resampled)

    #from sklearn.preprocessing import RobustScaler
    #scaler = RobustScaler()
    #X_resampled = scaler.fit_transform(X_resampled)

    #from sklearn.preprocessing import Normalizer
    #scaler = Normalizer()
    #X_resampled_smote = scaler.fit_transform(X_resampled_smote)


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
    #Train_X, Test_X, Train_y, Test_y = model_selection.train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=2020)

    import numpy as np
    from sklearn.linear_model import LogisticRegression, SGDClassifier
    from sklearn.svm import SVC
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score
    from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
    from xgboost import XGBClassifier
    import lightgbm as lgb #lgb.LGBMClassifier
    from catboost import CatBoostClassifier
    from sklearn.neighbors import KNeighborsClassifier


    """
    def classifier():
        scores=[]
        f1s=[]
        aucs=[]
        precisions=[]
        recalls=[]
        models = [
            LogisticRegression(class_weight='balanced',random_state=2020,max_iter=1000),
            SGDClassifier(class_weight='balanced',random_state=2020),
            SVC(class_weight='balanced',random_state=2020), #RBF
            DecisionTreeClassifier(class_weight='balanced',random_state=2020),
            RandomForestClassifier(class_weight='balanced',random_state=2020),
            ExtraTreesClassifier(class_weight='balanced',random_state=2020),
            #XGBClassifier(random_state=2020,n_estimators=200,eta=0.1, max_depth=6, objective='binary:logistic',scale_pos_weight = float(WT / AGO)),
            #lgb.LGBMClassifier(class_weight='balanced', random_state=2020),
            XGBClassifier(random_state=2020, n_estimators=200, eta=0.1, max_depth=6, objective='binary:logistic',scale_pos_weight=float(Weight)),
            lgb.LGBMClassifier(random_state=2020, n_estimators=100, max_depth=6,num_leaves=64, scale_pos_weight= float(Weight), objective='binary') #,class_weight='balanced'
            #CatBoostClassifier(class_weights={0: float(WT / All), 1: float(AGO / All)}, silent=True, random_state=2020,iterations=100, max_depth=None)
                ]

        for model in models:
            score = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='accuracy')
            mean_score = np.round(np.mean(score),5)
            scores.append(mean_score)

            precision = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='precision')
            mean_precision = np.round(np.mean(precision),5) #https://blog.csdn.net/liuweiyuxiang/article/details/100574386
            precisions.append(mean_precision)

            recall = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='recall')
            mean_recall = np.round(np.mean(recall),5) #https://blog.csdn.net/liuweiyuxiang/article/details/100574386
            recalls.append(mean_recall)

            f1 = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='f1')
            mean_f1 = np.round(np.mean(f1),5) #https://blog.csdn.net/liuweiyuxiang/article/details/100574386
            f1s.append(mean_f1)

            auc = cross_val_score(model, X=Train_X, y=Train_y, cv=10, scoring='roc_auc')
            mean_auc = np.round(np.mean(auc), 5)  # https://blog.csdn.net/liuweiyuxiang/article/details/100574386
            aucs.append(mean_auc)

        return scores, f1s, aucs, precisions, recalls


    scores, f1s, aucs, precisions, recalls = classifier()
    #print('accurcay', scores)
    #print('precision', precisions)
    #print('recall', recalls)
    #print('f1', f1s)
    #print('auc', aucs)



    print(str(scores).replace("[","").replace("]","")+"\t"+
          str(precisions).replace("[","").replace("]","")+"\t"+
          str(recalls).replace("[","").replace("]","")+"\t"+
          str(f1s).replace("[","").replace("]","")+"\t"+
          str(aucs).replace("[","").replace("]",""))

    """






    #Model
    lgb=XGBClassifier(random_state=2020,n_estimators=200,eta=0.1, max_depth=6, objective='binary:logistic',scale_pos_weight = float(Weight))  ##class_weight='balanced',is_unbalance=True，
    #lgb = lgb.LGBMClassifier(random_state=2020, n_estimators=100, max_depth=6,num_leaves=64,scale_pos_weight= float(Weight), objective='binary') #,class_weight='balanced'





    model=lgb.fit(Train_X, Train_y)


    print("Model Score -> ", model.score(Train_X, Train_y) * 100)

    # predict the labels on validation dataset
    predictions = model.predict(Test_X)

    # Use accuracy_score function to get the accuracy
    from sklearn.metrics import accuracy_score

    print("Accuracy Score -> ", accuracy_score(predictions, Test_y) * 100)

    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import classification_report

    confusion_matrix(predictions, Test_y)

    print(classification_report(predictions, Test_y))

    #模型生成
    import pickle

    f = open('model2.pickle', 'wb')
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
        train_sizes, train_scores, test_scores = learning_curve(
            estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
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

    plot_learning_curve(lgb, title, Train_X, Train_y, (0.8, 1.01), cv=cv, n_jobs=1)
    plt.savefig('lc2.png')
    #plt.show()