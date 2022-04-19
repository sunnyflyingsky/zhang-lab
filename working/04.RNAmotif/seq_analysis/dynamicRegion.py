from IPython.core.display import display, HTML
from numpy.lib import diff
from traitlets.traitlets import CInt
#display(HTML("<style>.container { width:90% !important; }</style>"))
import itertools
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binom
import pandas as pd
import GAP

import os
import sys
import getopt

def dynamicPipeline(inputDict:dict={},outputFolder:str='',refFolder:str='ref',\
    window1:int=101,step1:int=1,window2:int=30,step2:int=10,\
    cutoff:float=0.0001,flag:bool=True,start_step:int=5,c1:float=0.6,\
    c2:int=0,r1:str='',r2:str='',dyn:int=5,min_num:int=60):
    resPath1 = outputFolder+'/results.txt'
    resPath2 = outputFolder+'/results_filter.txt'
    resPath3 = outputFolder+'/results_merge_filter.txt'
    resPath4 = outputFolder+'/results_merge_filter_window.txt'
    resPath5 = outputFolder+'/results_merge_filter_window_anno.txt'
    resPath6 = outputFolder+'/results_merge_filter_window_anno_utr3.txt'
    resPath7 = outputFolder+'/window-anno.bed'
    
    print('#'*100,end='\n\n')
    if start_step<=1:
        NES_R1 = read(inputDict["NES"][0])
        NES_R2 = read(inputDict["NES"][1])
        ERM_R1 = read(inputDict["ERM"][0])
        ERM_R2 = read(inputDict["ERM"][1])
        NES_overlap = overlap(NES_R1, NES_R2)
        ERM_overlap = overlap(ERM_R1, ERM_R2)
        print('组内相交的转录组: ',end='')
        print(NES_overlap,ERM_overlap)
        In_diff = Step01(ERM_R1,ERM_R2,NES_R1,NES_R2)
        print('组内差异: ',end='')
        print(In_diff)
        del NES_R1,NES_R2,ERM_R1,ERM_R2,NES_overlap,ERM_overlap

        NES = read(inputDict["NES"][-1])
        ERM = read(inputDict["ERM"][-1])
        Step02(NES,ERM)
        Step03(NES,ERM,resPath1,window1,step1,In_diff,dyn,min_num)
        del NES,ERM
    
    print('#'*100,end='\n\n')
    if start_step<=2:
        avgPath = refFolder+'/avg.txt'
        avg_diff = Step04(resPath1,avgPath)
        #Step05(resPath1,avg_diff,resPath2,cutoff)
        print('平均信号0.9分位数: ',end='')
        print(avg_diff)
    
    print('#'*100,end='\n\n')
    if start_step<=3:
        Step06(resPath1,resPath3)
        #Step06(resPath2,resPath3)
    print('#'*100,end='\n\n')
    if start_step<=4:
        Step07(resPath3,resPath4,window2,step2)
    print('#'*100,end='\n\n')
    if start_step<=5:
        Step08(resPath4,resPath5,r1,r2,c1=c1,c2=c2)
    print('#'*100,end='\n\n')
    if start_step<=6:
        Step09(resPath5,resPath6,resPath7,flag)
    print('#'*100,end='\n\n')

def Step01(ERM_R1,ERM_R2,NES_R1,NES_R2):
    '''01 组内差异计算'''
    print('step01 start(10%):')
    ERM_dis,ERM_diff = calculate_distribution(ERM_R1,ERM_R2,c='b',figSavePath='')
    NES_dis,NES_diff = calculate_distribution(NES_R1,NES_R2,c='r',figSavePath='')
    Ingroup = NES_dis + ERM_dis
    In_diff = get_distribution(Ingroup,0.95,c='y',figSavePath='')
    #print(ERM_diff,NES_diff,In_diff)
    return In_diff

def Step02(NES,ERM):
    '''02 组间差异计算'''
    print('step02 start(20%):')
    Overgroup_overlap = overlap(NES,ERM)
    Over_dis,Over_diff = calculate_distribution(NES,ERM,c='g',figSavePath='')
    print('组间差异: ',end='')
    print(Overgroup_overlap,Over_diff)

def Step03(NES,ERM,resPath1,window1,step1,In_diff,dyn=5,min_num=60):
    '''03 滑动窗口扫描'''
    print('step03 start(30%):')
    print('dynamicWindows={}\ndynamicSteps={}\ndynamic_region_file=\"{}\"\n'.format(str(window1),str(step1),resPath1))
    fw = open(resPath1, 'w')
    title="id"+"\t"+"start"+"\t"+"end"+"\t"+"dyn_nul"+"\t"+"val_nul"+"\t"+"avg"+"\t"+"label"+"\t"+"p_value"+"\n"
    fw.write(title)
    ws=window1
    ss=step1
    for name in ERM.keys():
        if (name in NES.keys()):
            shape1 = NES[name]
            shape2 = ERM[name]
            trans, pos, dyn_nul, val_nul, avg, label, pvalue = \
                dynamic(name, shape1, shape2, window=ws, step=ss, dyn=dyn, min_num=min_num, score=In_diff)
            print('-',end='')
            for i in range(0,len(pos),ss):
                results = str(trans[i])+"\t"+str(pos[i])+"\t"+str(pos[i]+ws)+"\t"+ \
                    str(dyn_nul[i])+"\t"+str(val_nul[i])+"\t"+str(avg[i])+"\t"+str(label[i])+"\t"+str(pvalue[i])+"\n"
                fw.write(results)
    fw.close()
    print('|')

def Step04(resPath1,avgPath):
    '''04 获取avg signals的threshold'''
    print('step04 start(40%):')
    print('average_shape_score_file=\"{}\"\n'.format(avgPath))
    subprocess.call(["awk \'{print $6}\' %s > %s"%(resPath1,avgPath)],shell=True)
    with open(avgPath) as f:
        next(f) #挑过第一行
        lines = f.read().splitlines()
    lines = list(map(float, lines))
    avg_diff = get_distribution(lines,0.9,c='b',figSavePath='')
    return avg_diff
    
def Step05(resPath1,avg_diff,resPath2,cutoff):
    '''05 数据过滤'''
    print('step05 start(50%):')
    print('filtered_dynamic_region_file=\"{}\"\nfilter_cutoff={}\n'.format(resPath2,str(cutoff)))
    df=pd.read_table(resPath1,sep='\t')
    df1 = df[(df['avg'] > avg_diff) & (df['p_value'] < cutoff)]
    df1.to_csv(resPath2, sep='\t', index=False)
    
def Step06(resPath2,resPath3):
    '''06 对overlap的region进行merge操作'''
    print('step06 start(60%):')
    print('merged_filtered_dynamic_region_file=\"{}\"\n',format(resPath3))
    df=pd.read_table(resPath2,sep='\t')
    df_merge = merge_intervals(df, id_col="id", start_col="start", \
        stop_col="end", dyn_col="dyn_nul", val_col="val_nul", \
        avg_col="avg", label_col="label", pvalue_col="p_value")
    df_merge.to_csv(resPath3, sep='\t', index=False)
    
def Step07(resPath3,resPath4,window2,step2):
    '''07 对merge后的大region进行滑动窗口切割'''
    print('step07 start(70%):')
    print('slidingWindows={}\nslidingSteps={}\nslided_merged_filtered_dynamic_region_file=\"{}\"\n'.format(\
        str(window2),str(step2),resPath4))
    subprocess.call(["python slide_window.py -i %s -w %s -s %s > %s"%(resPath3,window2,step2,resPath4)],shell=True)

def Step08(resPath4,resPath5,r1,r2,c1,c2):
    '''08 注释文件'''
    print('step08 start(80%):')
    print('ref_transcriptome_file=\"{}\"\nref_genomeCoor_file=\"{}\"\nNone_cutoff={}\nstruct_cutoff={}\nanno_slided_merged_filtered_dynamic_region_file=\"{}\"\n'.format(\
        r1,r2,str(c1),str(c2),resPath5))
    df_merge=pd.read_table(resPath4,sep='\t')
    homo_parser = GAP.init(r2, r1)
    df_merge = annotation(df_merge,homo_parser,c1,c2)
    df_merge['label']
    df_merge.to_csv(resPath5, sep='\t', index=False)

def Step09(resPath5,resPath6,resPath7,flag):
    '''09 获取UTR3区域的dynamic region'''
    print('step09 start(90%):')
    if flag:
        print('UTR3_region_file=\"{}\"\nfinal_file=\"{}\"\nflag={}\n'.format(resPath6,resPath7,flag))
        subprocess.call(["awk \'$1==\"id\"{print $0}\' %s > %s"%(resPath5,resPath6)],shell=True)
        subprocess.call(["awk \'$9==\"utr3\"{print $0}\' %s >> %s"%(resPath5,resPath6)],shell=True)
    else:
        print('final_AllRegion_file=\"{}\"\nflag={}\n'.format(resPath7,flag))
        resPath6 = resPath5
    subprocess.call(["awk \'{print $1 \"\t\" $2 \"\t\" $3 \"\t\" $1 \"\t\" $12 \"\t\" $13 \"\t\" $11 \"\t\" $9 \"\t\" $7 \"\t\" $19 \"\t\" $20 \"\t\" $21 \"\t\" \"*\"}\' %s > %s"%(resPath6,resPath7)],shell=True)

def calculate_distribution(R1:dict,R2:dict,a:float=0.95,c=None,figSavePath:str=''):
    '''
    '''
    merge=[]
    for name in R1.keys():
        if (name in R2.keys()):
            shape1 = R1[name]
            shape2 = R2[name]
            merge.append(diff_abs(shape1, shape2))
    dis =list(itertools.chain.from_iterable(merge))
    diff = get_distribution(dis,a,c,figSavePath)
    return dis,diff

def get_distribution(dis:list,a:float=0.95,c=None,figSavePath:str=''):
    '''
    plot distribution fig, get the 100a% number of distribution
    '''
    diff = np.quantile(dis, a)
    #sns.histplot(dis,color=c,alpha=0.2)
    #plt.savefig(figSavePath)
    return diff

def read(shape_file:str):
    '''
    read shape file
    '''
    data = {}
    f = open(shape_file, 'r')
    i = 0
    for line in f:
        i = i + 1
        line = line.strip('\n')
        sent = line.split('\t')
        #sent1 = sent[0].split('.') 
        data[sent[0]] = sent[3:]
    plen = i
    return data
    f.close()

def overlap(shape1:dict, shape2:dict):
    '''
    get the overlap transcripts of two shape info dict
    '''
    num=0
    for name in shape1.keys():
        if (name in shape2.keys()):
            num += 1
    return num

def diff_abs(Shape1:dict, Shape2:dict):
    '''
    get the difference value of two shape info dict
    return abs(shape1-shape2)
    '''
    minus=[]
    if len(Shape1) != len(Shape2):
        pass
    else:
        for i in range(len(Shape1)):
            if (Shape1[i] != 'NULL') and (Shape2[i] != 'NULL'):
                minus.append(abs(float(Shape1[i]) - float(Shape2[i])))
        return minus

def diff_(Shape1:dict, Shape2:dict):
    '''
    get the difference value of two shape info dict
    return shape1-shape2
    if shape1>shape2, you will get positive value,
    if shape1<shape2, you will get negative value,
    '''
    minus=[]
    if len(Shape1) != len(Shape2):
        pass
    else:
        for i in range(len(Shape1)):
            if (Shape1[i] != 'NULL') and (Shape2[i] != 'NULL'):
                minus.append(float(Shape1[i]) - float(Shape2[i]))
        return minus

def dynamic(Name, Shape1:dict, Shape2:dict, window:int=30, step:int=1,\
    score:float=0.5, dyn:int=0, min_num:int=20):
    '''
    calculate dynamic region with a sliding windows
    '''
    test_est = []
    test_p = []
    ll_dist1 = []
    pos=[]
    pvalue=[]
    trans=[]
    label = []
    j=0
    if len(Shape1) != len(Shape2):
        pass
    else:
        for i in range(0,len(Shape1) - window,step):
            a = Shape1[i:(i+window)]
            b = Shape2[i:(i+window)]
            minus = []
            dyn_nul = 0
            val_nul = 0
            seq = ''
            for it in range(len(a)):
                if (a[it] != 'NULL') and (b[it] != 'NULL'):
                    minus.append(abs(float(a[it]) - float(b[it])))
                    if (abs(float(a[it]) - float(b[it])) >= score):
                        dyn_nul = dyn_nul + 1
                        if float(a[it])>float(b[it]):
                            seq+='+'
                        else:
                            seq+='-'
                    else:
                        seq+='*'
                    val_nul = val_nul + 1
                else:
                    seq+='.'
            j+=1
            x = np.array(minus)
            if val_nul >= min_num:
                minus_avg = round(x.mean(),5)
            else:
                dyn_nul = -1
                val_nul = -1
                minus_avg = -1
              
            if dyn_nul>=dyn:
                trans.append(Name)
                pos.append(j)
                test_est.append(dyn_nul)
                test_p.append(val_nul)
                ll_dist1.append(minus_avg)
                pvalue.append(round(1-binom.cdf(dyn_nul, val_nul, 0.05),5))
                label.append(seq)
            else:
                pass
   
    return trans, pos, test_est, test_p, ll_dist1, label, pvalue

def merge_intervals(df, id_col="id", start_col="start", \
    stop_col="end", dyn_col="dyn_nul", val_col="val_nul", \
    avg_col="avg", label_col="label", pvalue_col="p_value"):
    '''
    merge dynamic region
    '''
    out = []
    #df = df.sort_values(start_col)
    first_row = True
    for ix, row in df.iterrows():
        if first_row:
            current_start = row[start_col]
            current_stop = row[stop_col]
            current_id = row[id_col]
            current_dyn = row[dyn_col]
            current_val = row[val_col]
            current_avg = row[avg_col]
            current_label = row[label_col]
            current_pvalue = row[pvalue_col]
            first_row = False
        start = row[start_col]
        stop = row[stop_col]
        id = row[id_col]
        dyn = row[dyn_col]
        val = row[val_col]
        avg = row[avg_col]
        label = row[label_col]
        pvalue = row[pvalue_col]
        if id == current_id:
            if start > stop:
                raise Exception("Starts must be greater stops!")
            if start > current_stop:
                # Gap between segments: output current segment and start a new one.
                out.append((current_id, current_start, current_stop, current_dyn, current_val, current_avg, current_label, current_pvalue))
                current_id, current_start, current_stop,current_dyn,current_val, current_avg, current_label, current_pvalue = \
                    id, start, stop, dyn, val, avg, label, pvalue
            else:
                # Segments adjacent or overlapping: merge.
                if current_stop>=stop:
                    current_label = current_label[:]
                else:
                    current_label = current_label[:]+label[current_stop-start:]
                current_stop = max(current_stop, stop)
                #current_dyn = max(current_dyn, dyn) #获得对应的索引
        elif id != current_id:
            if start > stop:
                raise Exception("Starts must be greater stops!")
            #if start <= current_stop:
                # Gap between segments: output current segment and start a new one.
            out.append((current_id, current_start, current_stop, current_dyn, current_val, current_avg, current_label, current_pvalue))
            current_id, current_start, current_stop,current_dyn,current_val, current_avg, current_label, current_pvalue = \
                id, start, stop, dyn, val, avg, label, pvalue
            #else:
                #print(id)
                # Segments adjacent or overlapping: merge.
                #current_stop = max(current_stop, stop)
                #if current_stop>=stop:
                #    current_label = current_label[:]
                #else:
                #    current_label = current_label[:]+label[current_stop-start:]
    out.append((current_id, current_start, current_stop, current_dyn, current_val, current_avg, current_label, current_pvalue))
    return pd.DataFrame(out, columns=[id_col, start_col, stop_col, dyn_col, val_col, avg_col, label_col, pvalue_col])

def annotation(df_merge:dict,homo_parser,c1:float=0.6,c2:int=0):
    '''
    annotation
    '''
    region=[]
    chrom=[]
    strand=[]
    pos_s=[]
    pos_e=[]
    name=[]
    identity=[]
    category=[]
    length=[]
    seq=[]
    flags=[]
    less = []
    more = []
    for index, row in df_merge.iterrows():
        #获得转录本位置信息
        ft = homo_parser.getTransFeature(row['id'])
        utr5_s, utr5_e, cds_s, cds_e, utr3_s, utr3_e = ft['utr_5_start'],ft['utr_5_end'], ft['cds_start'], ft['cds_end'], ft['utr_3_start'],ft['utr_3_end']
        if (row['end']<=utr5_e): #row['start']>=utr5_s and 
            region.append("utr5")
        elif(row['start']>=cds_s and row['end']<=cds_e):
            region.append("cds")
        elif(row['start']>=utr3_s and row['end']<=utr3_e): 
            region.append("utr3")
        elif(row['start']>=utr5_s and row['start']<=utr5_e and row['end']>=cds_s and row['end']<=cds_e):
            region.append("utr5+cds")
        elif(row['start']>=cds_s and row['start']<=cds_e and row['end']>=utr3_s and row['end']<=utr3_e):
            region.append("cds+utr3")
        else:
            region.append("unknown")
        label = row['label']
        if label.count('.')/len(label)>=c1:
            flags.append('unknown')
            less.append(0)
            more.append(0)
        else:
            more.append(label.count('+'))
            less.append(label.count('-'))
            if label.count('+')-label.count('-')>c2:
                flags.append('more')
            elif label.count('+')-label.count('-')<-c2:
                flags.append('less')
            else:
                flags.append('unknown')
        chr, std, s, e, gene_name, gene_id, gene_type, trans_len = ft['chr'],ft['strand'],ft['start'],ft['end'],ft['gene_name'],ft['gene_id'], ft['gene_type'], ft['trans_len']
            
        chrom.append(chr)
        strand.append(std)
        pos_s.append(s)
        pos_e.append(e)
        name.append(gene_name)
        identity.append(gene_id)
        category.append(gene_type)    
        length.append(trans_len)
        
        GAPDH = homo_parser.getTransSeq(row['id'])
        seq.append(GAPDH[row['start']:row['end']])

    df_merge['region'] = region
    df_merge['chr'] = chrom
    df_merge['strand'] = strand
    df_merge['gent_start'] = pos_s
    df_merge['gene_end'] = pos_e
    df_merge['gene_name'] = name
    df_merge['gene_id'] = identity
    df_merge['gene_type'] = category
    df_merge['trans_len'] = length
    df_merge['seq'] = seq
    df_merge['flag'] = flags 
    df_merge['more'] = more
    df_merge['less'] = less

    return df_merge

#Usage = ''
Usage = 'Usage: ' + sys.argv[0]
Usage += '''
    <Requires>
    -i parameter.in

    you should input the input information in parameter.in, \n
    including: 
    inputfolder, all of your input file should be in this folder path, including eg1,eg2,egF,cg1,cg2,cgF.
    refFolder, all of your annotation file should be in this folder path, including genomeCoor and transcriptome file.
    ref_genomeCoor_file, genomeCoor file.
    ref_transcriptome_file, transcriptome file.
    outputFolder, we will output the results to this path. and if you start from step>1(if you first use, please set the start_step to 1), please put the file in this folder, including dynamic_region_demo, filtered, merged, anno, etc.  
    inputEgFile1, experimental icshape score, group1, for calculate difference.
    inputEgFile2, experimental icshape score, group2, for calculate difference.
    inputCgFile1, control icshape score, group1, for calculate difference.
    inputCgFile2, control icshape score, group2, for calculate difference.
    inputEgFileFinal, experimental icshape score for calculate dynamic region
    inputCgFileFinal, control icshape score for calculate dynamic region

    [Options]
    dynamicWindows, window size for dynamic search
    dynamicSteps, window step for dynamic search
    slidingWindows, slide window size for cutting the merged region to get same size regions 
    slidingSteps, slide window step for cutting the merged region to get same size regions
    filter_cutoff, filter cutoff for P-value
    None_cutoff, cutoff for label the region with much None score siganl. the larger it is, the more region will be labeled unknown.  
    struct_cutoff, cutoff for label the region to more single-structure or less single-structure. the larger it is, the more region will be labeled unknown. 
    start_tep, which step you want to start. if you first use, please set it to 1. it should be a int number from 1 to 6.  
                1 means from shape file to get all the results file
                2 means from dynamic_region_demo file to get onther results file
                3 means from dynamic_region_filtered file to get following results file
                4 means from dynamic_region_filtered_merged file to get following results file
                5 means from dynamic_region_filtered_merged_slided file to get following results file
                6 means from dynamic_region_filtered_merged_slided_anno file to get final window-anno.bed file
    dynamic_positive_num, cutoff for select dynamic windows, which have more sites(score>indiff) than this num. 
    dynamic_valid_num, cutoff for select dynamic windows, which have more valid score than this num.
    UTR3_or_ALLREGION, if you want to get UTR3 region. True for UTR3, False for all region.
'''
Usage += "EX: python " + sys.argv[0] + ''' -i parameter.in'''
if len(sys.argv)<2 or not sys.argv[1].startswith('-'):sys.exit(Usage)

if __name__ == '__main__':
    inpath = 'parameter.in'

    oplist,alist = getopt.getopt(sys.argv[1:],'hi:')
    for opt in oplist:
        if opt[0] == '-h':sys.exit(Usage)
        elif opt[0] == '-i':inpath = opt[1]
        else: 
            sys.exit(Usage)

    inputFolder = 'data'
    outputFolder = 'res'
    refFolder = 'ref'
    window1=101
    step1=1
    window2=30
    step2=10
    cutoff=0.0001
    start_step=1
    flag=True
    c1=0.6
    c2=0
    dyn=5
    min_num=60
    r1='ref/hg38_transcriptome.fa'
    r2='ref/hg38.genomeCoor.bed'

    inputDict = {'NES':['','',''],\
        'ERM':['','','']}
    
    #获取对应的参数信息
    with open(inpath) as read_obejct:
        for line in read_obejct:
            if line.startswith('inputFolder'):
                inputFolder=line.strip().split(' ')[-1]
            elif line.startswith('refFolder'):
                refFolder=line.strip().split(' ')[-1]
            elif line.startswith('outputFolder'):
                outputFolder=line.strip().split(' ')[-1]
            elif line.startswith('inputEgFile1'):
                eg1 = line.strip().split(' ')[-1]
                inputDict['NES'][0]=inputFolder+'/'+eg1
            elif line.startswith('inputEgFile2'):
                eg2 = line.strip().split(' ')[-1]
                inputDict['NES'][1]=inputFolder+'/'+eg2
            elif line.startswith('inputCgFile1'):
                cg1 = line.strip().split(' ')[-1]
                inputDict['ERM'][0]=inputFolder+'/'+cg1
            elif line.startswith('inputCgFile2'):
                cg2 = line.strip().split(' ')[-1]
                inputDict['ERM'][1]=inputFolder+'/'+cg2
            elif line.startswith('inputEgFileFinal'):
                egf = line.strip().split(' ')[-1]
                inputDict['NES'][-1]=inputFolder+'/'+egf
            elif line.startswith('inputCgFileFinal'):
                cgf= line.strip().split(' ')[-1]
                inputDict['ERM'][-1]=inputFolder+'/'+cgf
            elif line.startswith('ref_genomeCoor_file'):
                r2 = refFolder+'/'+line.strip().split(' ')[-1]
            elif line.startswith('ref_transcriptome_file'):
                r1 = refFolder+'/'+line.strip().split(' ')[-1]
            elif line.startswith('dynamicWindows'):
                try:
                    window1=abs(int(line.strip().split(' ')[-1]))
                except:
                    print('dynamicWindows should be int!')
                    window1=101
            elif line.startswith('dynamicSteps'):
                try:
                    step1=abs(int(line.strip().split(' ')[-1]))
                except:
                    print('dynamicSteps should be int!')
                    step1=1
            elif line.startswith('slidingWindows'):
                try:
                    window2=abs(int(line.strip().split(' ')[-1]))
                except:
                    print('slidingWindows should be int!')
                    window2=30
            elif line.startswith('slidingSteps'):
                try:
                    step2=abs(int(line.strip().split(' ')[-1]))
                except:
                    print('slidingSteps should be int!')
                    step2=10
            elif line.startswith('filter_cutoff'):
                try:
                    cutoff=float(line.strip().split(' ')[-1])
                except:
                    print('filter_cutoff should be float!')
                    cutoff=0.0001
            elif line.startswith('None_cutoff'):
                try:
                    c1=float(line.strip().split(' ')[-1])
                except:
                    print('None_cutoff should be float!')
                    c1=0.6
            elif line.startswith('struct_cutoff'):
                try:
                    c2=abs(int(line.strip().split(' ')[-1]))
                except:
                    print('struct_cutoff should be int!')
                    c2=0
            elif line.startswith('start_step'):
                try:
                    start_step=int(line.strip().split(' ')[-1])
                except:
                    print('start_step should be int!')
                    start_step=1
                if start_step<1 or start_step>=7:
                    print('start_step should be int number in range[1,6]')
            elif line.startswith('dynamic_positive_num'):
                try:
                    dyn=abs(int(line.strip().split(' ')[-1]))
                except:
                    print('dynamic_positive_num should be int!')
                    dyn=5
            elif line.startswith('dynamic_valid_num'):
                try:
                    min_num=abs(int(line.strip().split(' ')[-1]))
                except:
                    print('dynamic_valid_num should be int!')
                    min_num=60
            elif line.startswith('UTR3_or_ALLREGION'):
                flag=line.strip().split(' ')[-1]
                if flag=='True' or flag=='1' or flag=='true' or flag=='TRUE':
                    flag=True
                else:
                    flag=False

    for key,value in inputDict.items():
        for p in value:
            try:
                open(p)
            except:
                print('filepath "'+p+'" is not exist!')
                sys.exit(Usage)
    print('run(0%)!')
    dynamicPipeline(inputDict=inputDict,outputFolder=outputFolder,refFolder=refFolder,\
        window1=window1,step1=step1,window2=window2,step2=step2,cutoff=cutoff,\
        flag=True,start_step=start_step,c1=c1,c2=c2,\
        r1=r1,r2=r2,dyn=dyn,min_num=min_num)
    print('finished(100%)!')
