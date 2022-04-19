#import itertools
import subprocess
import numpy as np
import pandas as pd
import sys
import os
import getopt
from threading import Thread
from time import sleep, ctime
#import GAP
#import General
#hg38_parser = GAP.init("/data/zhangjiasheng/data/ref_info/hg38.genomeCoor.bed", "/data/zhangjiasheng/data/ref_info/hg38_transcriptome.fa")

def sliding(seq,w=200,s=200):
    '''
    '''
    seq_list = []
    for i in range(0,len(seq)-w,s):
        seq_list.append(seq[i:i+w])
    return seq_list

def shape2structure(inputpath,outputpath,index=1):
    seq_list = {}
    with open(inputpath) as read_object:
        i = 0
        for line in read_object:
            i+=1
            info = line.strip().split('\t')
            if index == 1:
                seq_list[info[0]]=sliding(info[3:],200,200)
            elif index==0:
                pass
                #seq_list[info[0]]=sliding(hg38_parser.getTransSeq(info[0]),200,200)
    with open(outputpath,'w') as write_object:
        for key,value in seq_list.items():
            for i,seq in enumerate(value):
                write_object.write('>{}\n'.format(key+'_'+str(i)))
                write_object.write('{}\n'.format(seq))

def one_compute_by_one(inpath,refpath):
    '''
    从分别制作为fasta格式的seq和shape数据中计算结构
    获取200bp序列以及它的点括号图
    '''
    shape = {}
    with open(refpath) as read_object:
        for line in read_object:
            if line.startswith('>'):
                name = line.strip().replace('>','')
            else:
                try:
                    shape[name].append(line.strip())
                except:
                    shape[name] = []
    seq = {}
    with open(inpath) as read_object:
        for line in read_object:
            if line.startswith('>'):
                name=line.strip().replace('>','')
            else:
                try:
                    seq[name].append(line.strip())
                except:
                    seq[name] = []
    for key,value in seq.items():
        for i in range(len(value)):
            fapath = '/data/zhangjiasheng/data/shape2structure/01_shape/tmp/fa.fasta'
            oupath = '/data/zhangjiasheng/data/shape2structure/01_shape/tmp/temp.ct'
            shpath = '/data/zhangjiasheng/data/shape2structure/01_shape/tmp/sh.shape'
            with open(fapath,'w') as write_object:
                write_object.write('>'+key+'\n')
                write_object.write(value[i])
            with open(shpath,'w') as write_object:
                shape_number = shape[key][i].split(',')
                for j,s in enumerate(shape_number):
                    if s.strip()=='NULL':
                        write_object.write('{} {}\n'.format(str(j+1),'-999'))
                    else:
                        write_object.write('{} {}\n'.format(str(j+1),s))
            subprocess.call(["Fold %s %s -sh %s"%(fapath,oupath,shpath)],shell=True)
            subprocess.call(["ct2dot %s %s %s"%(oupath,1,oupath.replace('temp.ct','temp.dot'))],shell=True)
            subprocess.call(["cat %s >> /data/zhangjiasheng/data/shape2structure/02_structure/%s.dot"%(oupath.replace('temp.ct','temp.dot'),'ERM')],shell=True)

def concat(inpath,outpath):
    '''
    '''
    name = ''
    structure_data = {}
    with open(inpath) as read_object:
        for line in read_object:
            name_current = line.strip().split(' ')[-1]
            if line.startswith('>'):
                if name_current!= name:
                    name = name_current
                    structure_data[name] = ['','']
                else:
                    pass
            elif 'A' in line or 'T' in line or 'C' in line or 'G' in line:
                structure_data[name][0]+=line.strip()
            else:
                structure_data[name][1]+=line.strip()
    with open(outpath,'w') as write_object:
        for name,info in structure_data.items():
            write_object.write('>{}\n'.format(name))
            write_object.write('{}\n'.format(info[0]))
            write_object.write('{}\n'.format(info[1]))

def motif_search(inpath,tmppath,outpath):
    '''
    '''
    threadList = []
    structure_data = {}
    i = 0
    with open(inpath) as read_object:
        for line in read_object:
            if line.startswith('>'):
                name = line.strip().replace('>','')
                structure_data[name]=['','']
            elif 'A' in line or 'T' in line or 'C' in line or 'G' in line:
                structure_data[name][0]+=line.strip()
            else:
                structure_data[name][1]+=line.strip()
                if len(structure_data.keys())>=200:
                    i+=1
                    subprocess.call(["mkdir /data/zhangjiasheng/data/shape2structure/03_single_fa_dot/tmp"+str(i)],shell=True)
                    threadList.append(Thread(target=motif_search_one_thread,\
                        args=(structure_data.copy(),tmppath+'tmp'+str(i)+'/',outpath.replace('.txt',str(i)+'.txt'),i)))
                    structure_data = {}
    for t in threadList:
        t.start()
    for t in threadList:
        t.join()
                            
def motif_search_one_thread(structure_data,tmppath,outpath,t_index):
    f = open(outpath,'w') 
    f.close()
    filelist = ['single_dot_final.out','single_dot_mfinal.out','single_dot_rfinal.out']
    for key,value in structure_data.items():
        for j in range(0,len(value[0])-31,15):
            name = key+'_'+str(j)
            seq_single = value[0][j:j+31]
            struct_single = value[1][j:j+31]
            f = open(tmppath+'single_dot.txt','w')
            f.write('>{}\n'.format(name))
            f.write('{}\n'.format(seq_single))
            f.write('{}\n'.format(struct_single))
            f.close()
            subprocess.call(["cd /data/zhangjiasheng/toolkits/RNALigands && perl run.pl -s %s"%(tmppath+'single_dot.txt')],shell=True)
            subprocess.call(["cd /data/zhangjiasheng/data/shape2structure/03_single_fa_dot/tmp"+str(t_index)+" && rm *4.out *3.out *2.out *1.out"],shell=True)
            for file_out in filelist:
                with open(tmppath+file_out) as read_object:
                    for line in read_object:
                        info = line.strip().split(';')
                        try:
                            drug = info[-3]
                            score = info[-1]
                            with open(outpath,'a+') as write_object:
                                write_object.write('{}\t{}\t{}\t{}\t{}\n'.format(name,drug,seq_single,struct_single,score))
                        except:pass

def data_filtered(inpath,refpath,outpath):
    '''
    '''
    drugs = {}
    with open(refpath) as read_object:
        for line in read_object:
            info = line.strip().split('\t')
            drugs[info[0]]=info[1]
    with open(outpath,'w') as write_object:
        with open(inpath) as read_object:
            for line in read_object:
                info = line.strip().split('\t')
                if info[1]=='3':
                    info[1]='Benzimidazole_3'
                try:
                    write_object.write('{}\t{}\t{}\t{}\n'.format(drugs[info[1]],info[2],info[3],info[4]))
                except:
                    print(info[1])


def main(inputpath,outputpath,shapepath):
    #hape2structure(inputpath,outputpath,0)
    #shape2structure(inputpath,shapepath,1)
    inpath = outputpath
    refpath = shapepath
    #one_compute_by_one(inpath,refpath)
    inpath = '/data/zhangjiasheng/data/shape2structure/02_structure/ERM.dot'
    outpath = '/data/zhangjiasheng/data/shape2structure/02_structure/ERM_concat.dot'
    #concat(inpath,outpath)
    inpath = outpath
    tmppath = '/data/zhangjiasheng/data/shape2structure/03_single_fa_dot/' 
    outpath = '/data/zhangjiasheng/data/shape2structure/04_RNA_drug/RNA_Ligand.txt'
    motif_search(inpath,tmppath,outpath)
    inpath=outpath
    refpath = '/data/zhangjiasheng/data/shape2structure/04_RNA_drug/drug2smile_corrected.txt'
    outpath = '/data/zhangjiasheng/data/shape2structure/04_RNA_drug/RNA_Ligand_valid_bak.txt'
    #data_filtered(inpath,refpath,outpath)


if __name__=='__main__':
    inputpath = '/data/zhangjiasheng/data/shape2structure/01_shape/Human_293T_smartSHAPE.out'
    outputpath = '/data/zhangjiasheng/data/shape2structure/01_shape/in.fasta'
    shapepath = '/data/zhangjiasheng/data/shape2structure/01_shape/in.shape'
    main(inputpath,outputpath,shapepath)






