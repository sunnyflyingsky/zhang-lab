
    from xgboost import XGBClassifier
    #from sklearn.svm import SVC



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
    #print('auc', aucs)






    #模型生成
    import pickle

    xgb = XGBClassifier(random_state=2020, n_estimators=200, eta=0.1, max_depth=6, objective='binary:logistic',scale_pos_weight=float(Weight))
    model=xgb.fit(Train_X, Train_y)
    f = open('model.pickle', 'wb')
    pickle.dump(model, f)
    f.close()
