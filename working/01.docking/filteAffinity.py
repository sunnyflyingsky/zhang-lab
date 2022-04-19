import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
from numpy import polyfit, poly1d
import pandas as pd

import os
import re
import math


def distribution(inputfolder='docking/filter/affinity_data',flag=1,feature='weight'):
    '''
    计算小分子weight核密度曲线分布
    1: initial
    2: filtered
    '''
    column_name = ["ID1","ID2","smile","seq","struct","site","frag","Affinity","ID3","Weight"]
    autodock_native_path = inputfolder + '/Native_results/AutoDock_results.txt'
    dock6_native_path = inputfolder + '/Native_results/dock6_results.txt'
    if flag==1:
        autodock_filter_path = inputfolder + '/Docking_results/extension_autodock_final_revised.txt'
        dock6_filter_path = inputfolder + '/Docking_results/extension_dock6_final_revised.txt'
    elif flag==2:
        autodock_filter_path = 'docking/filter/res/autodock.csv'
        dock6_filter_path = 'docking/filter/res/dock6.csv'
    else:
        return
    
    native_data = pd.read_csv(autodock_native_path,sep='\t')
    filter_data = pd.read_csv(autodock_filter_path,sep='\t',names=column_name)
    native_data = native_data[native_data['Weight']>=200]
    native_data = native_data[native_data['Weight']<=600]
    filter_data = filter_data[filter_data['Weight']<=600]
    print('native: ',min(native_data['Weight']),max(native_data['Weight']))
    print('filter: ',min(filter_data['Weight']),max(filter_data['Weight']))
    #native_data.weight.plot(kind = 'hist', bins = 10, color = 'steelblue', edgecolor = 'black', density = True, label = '直方图')
    plotkdemap(native_data,filter_data,feature,'docking/filter/res/autodock_kde.png')

    native_data = pd.read_csv(dock6_native_path,sep='\t')
    filter_data = pd.read_csv(dock6_filter_path,sep='\t',names=column_name)
    native_data = native_data[native_data['Weight']>=200]
    native_data = native_data[native_data['Weight']<=600]
    filter_data = filter_data[filter_data['Weight']<=600]
    print('native: ',min(native_data['Weight']),max(native_data['Weight']))
    print('filter: ',min(filter_data['Weight']),max(filter_data['Weight']))
    plotkdemap(native_data,filter_data,feature,'docking/filter/res/dock6_kde.png')

def plotkdemap(native_data,filter_data,feature,respath):
    '''
    绘制核密度曲线
    '''
    plt.figure(figsize=(12,8))
    plt.rcParams.update({"font.size":20,'figure.dpi':300})
    plt.title("Kernel Density",fontsize=24)
    if feature=='weight':
        native_data.Weight.plot(kind = 'kde', color = 'red', label = 'native')
        filter_data.Weight.plot(kind = 'kde', color = 'blue', label = 'filter')
        plt.xlabel('molecular Weight')
    elif feature=='affinity':
        native_data.Affinity.plot(kind = 'kde', color = 'red', label = 'native')
        filter_data.Affinity.plot(kind = 'kde', color = 'blue', label = 'filter')
        plt.xlabel('molecular Affinity')
    else:
        return
    plt.ylabel('freq')
    plt.legend()
    plt.savefig(respath)
    plt.close()

def relationship(inputfolder='docking/filter/affinity_data',flag=1):
    '''
    绘制关系散点图
    '''
    column_name = ["ID1","ID2","smile","seq","struct","site","frag","Affinity","ID3","Weight"]
    autodock_native_path = inputfolder + '/Native_results/AutoDock_results.txt'
    dock6_native_path = inputfolder + '/Native_results/dock6_results.txt'
    if flag==1:
        autodock_filter_path = inputfolder + '/Docking_results/extension_autodock_final_revised.txt'
        dock6_filter_path = inputfolder + '/Docking_results/extension_dock6_final_revised.txt'
    elif flag==2:
        autodock_filter_path = 'docking/filter/res/autodock.csv'
        dock6_filter_path = 'docking/filter/res/dock6.csv'

    native_data = pd.read_csv(autodock_native_path,sep='\t')
    filter_data = pd.read_csv(autodock_filter_path,sep='\t',names=column_name)
    plotScatter(native_data, filter_data, respath='docking/filter/res/autodock_relationship.png')
    native_data = pd.read_csv(dock6_native_path,sep='\t')
    filter_data = pd.read_csv(dock6_filter_path,sep='\t',names=column_name)
    plotScatter(native_data, filter_data, respath='docking/filter/res/dock6_relationship.png')

def plotScatter(native_data,filter_data,respath):
    '''
    绘制关系散点图
    '''
    native_data = native_data[native_data['Affinity']<0]
    filter_data = filter_data[filter_data['Affinity']<0]
    plt.figure(figsize=(12,8))
    plt.rcParams.update({"font.size":20,'figure.dpi':300})
    plt.title("relationship",fontsize=24)
    plt.scatter(filter_data['Weight'],filter_data['Affinity'],color='blue',s=1,label='filter')
    plt.scatter(native_data['Weight'],native_data['Affinity'],color='red',s=1,label='native')
    plt.xlabel('molecular Weight')
    plt.ylabel('molecular Affinity')
    plt.legend()
    plt.savefig(respath)
    plt.close()

    plt.figure(figsize=(12,8))
    plt.rcParams.update({"font.size":20,'figure.dpi':300})
    plt.title("relationship",fontsize=24)
    plt.scatter(native_data['Weight'],native_data['Affinity'],color='red',s=3,label='native')
    plt.xlabel('molecular Weight')
    plt.ylabel('molecular Affinity')
    plt.legend()
    plt.savefig(respath.replace('.png','_native.png'))
    plt.close()

    plt.figure(figsize=(12,8))
    plt.rcParams.update({"font.size":20,'figure.dpi':300})
    plt.title("relationship",fontsize=24)
    plt.scatter(filter_data['Weight'],filter_data['Affinity'],color='blue',s=1,label='filter')
    plt.xlabel('molecular Weight')
    plt.ylabel('molecular Affinity')
    plt.legend()
    plt.savefig(respath.replace('.png','_filter.png'))
    plt.close()

def filter(inputfolder='docking/filter/affinity_data',method='bins'):
    '''
    过滤方法
    method: bins, flexible, linear, clusters,
    '''
    column_name = ["ID1","ID2","smile","seq","struct","site","frag","Affinity","ID3","Weight"]
    autodock_native_path = inputfolder + '/Native_results/AutoDock_results.txt'
    dock6_native_path = inputfolder + '/Native_results/dock6_results.txt'
    autodock_filter_path = inputfolder + '/Docking_results/extension_autodock_final_revised.txt'
    dock6_filter_path = inputfolder + '/Docking_results/extension_dock6_final_revised.txt'
    if method=='bins':
        native_data = pd.read_csv(autodock_native_path,sep='\t')
        filter_data = pd.read_csv(autodock_filter_path,sep='\t',names=column_name)
        filterByBins(native_data,filter_data,column_name,"docking/filter/res/autodock.csv")
        native_data = pd.read_csv(dock6_native_path,sep='\t')
        filter_data = pd.read_csv(dock6_filter_path,sep='\t',names=column_name)
        filterByBins(native_data,filter_data,column_name,"docking/filter/res/dock6.csv")
    elif method=='linear':
        native_data = pd.read_csv(autodock_native_path,sep='\t')
        filter_data = pd.read_csv(autodock_filter_path,sep='\t',names=column_name)
        filterByLinear(native_data,filter_data,column_name,'docking/filter/res/autodock.csv')
        native_data = pd.read_csv(dock6_native_path,sep='\t')
        filter_data = pd.read_csv(dock6_filter_path,sep='\t',names=column_name)
        filterByLinear(native_data,filter_data,column_name,'docking/filter/res/dock6.csv')
    elif method=='clusters':
        pass 

def filterByBins(native_data,filter_data,column_name,respath,cutoff=8):
    '''
    通过分bin的方式来过滤docking数据
    '''
    scale_Weight = [200,250,300,320,340,360,380,400,410,420,430,440,450,460,470,480,490,500,520,530,560,580,600,1500]
    print(len(filter_data['Affinity']))
    filter_data_scaled_all = pd.DataFrame(columns=column_name)
    for i in range(len(scale_Weight)-1):
        min_s,max_s=scale_Weight[i],scale_Weight[i+1]
        native_data_scaled = native_data[native_data['Weight']>=min_s]
        native_data_scaled = native_data_scaled[native_data_scaled['Weight']<max_s]
        native_data_scaled_l = native_data_scaled[native_data_scaled['RMSD']>=cutoff]
        if i<=340:p=0
        else:p=10
        try:threshold = np.percentile(native_data_scaled_l['Affinity'],p)
        except:continue
        native_data_scaled_r = native_data_scaled[native_data_scaled['RMSD']<=cutoff]
        native_data_scaled_r = native_data_scaled_r[native_data_scaled_r['Affinity']<threshold]
        try:threshold = np.percentile(native_data_scaled_r['Affinity'],100)
        except:continue

        print(threshold)
        filter_data_scaled = filter_data[filter_data['Weight']>=min_s]
        filter_data_scaled = filter_data_scaled[filter_data_scaled['Weight']<max_s]
        print(len(filter_data_scaled[filter_data_scaled['Affinity']<threshold]['Affinity']))
        filter_data_scaled = filter_data_scaled[filter_data_scaled['Affinity']<threshold]
        filter_data_scaled_all = filter_data_scaled_all.append(filter_data_scaled)
        #filter_data_scaled.to_csv("docking/filter/res/dock6_"+str(min_s)+"_"+str(max_s)+".csv",index=False,header=False,sep='\t')
    filter_data_scaled_all.to_csv(respath,index=False,header=False,sep='\t')

def filterByLinear(native_data,filter_data,column_name,respath):
    '''
    '''
    #print(native_data)
    param = polyfit(native_data['Weight'],native_data['Affinity'],1)
    delta_a = []
    for index,dataI in native_data.iterrows():
        affinity_p = linearFunc(dataI['Weight'],param[0],param[1])
        delta_a.append(dataI['Affinity'] - affinity_p)
    native_data['delta'] = delta_a
    native_data = native_data.sort_values(by='delta')
    L = 0
    l = 0
    for index,dataI in native_data.iterrows():
        if dataI['RMSD']<=8:
            print(index)
            param[1] = dataI['Affinity']-param[0]*dataI['Weight']
            L+=1
        elif l/L<=0.25:
            l+=1
        else:
            break
    filter_data_scaled_all = pd.DataFrame(columns=column_name)
    for index,dataI in filter_data.iterrows():
        if dataI['Affinity']<=param[0]*dataI['Weight']+param[1]:
            filter_data_scaled_all = filter_data_scaled_all.append(dataI)
            print(index)
    filter_data_scaled_all.to_csv(respath,index=False,header=False,sep='\t')

def filterByClusters(native_data,filter_data,column_name,respath):
    '''
    '''
    

def linearFunc(x,a,b):
    return a*x+b

if __name__ == '__main__':
    print('run!')
    filter(method='linear')
    distribution(flag=2,feature='weight')
    relationship(flag=2)
    
    














