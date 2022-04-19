import sys
import GAP
import General
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft,ifft,fftfreq
from matplotlib.pylab import mpl
import scipy.signal as signal

hg38_bed = "/Share2/home/zhangqf7/yuhan/0_annotation/Gencode/hg38.genomeCoor.bed"
hg38_seq = "/Share2/home/zhangqf7/yuhan/0_annotation/Gencode/hg38_transcriptome.fa"
hg38_parser = GAP.init(hg38_bed, hg38_seq)

def readFilterShape(inFile):
    SHAPE = General.load_shape(inFile)
    print("raw_size:", len(SHAPE))
    shape_keys = list(SHAPE.keys())
    for tid in shape_keys:
        length = len(SHAPE[tid])
        valid_shape = [ it for it in SHAPE[tid] if it != 'NULL' ]
        if 1.0*len(valid_shape)/length <= 0.2:
            del SHAPE[tid]
    
    print("new_size:", len(SHAPE))
    return SHAPE

IN_PATH = "./%s_FPKM.shape"

smart_1 = readFilterShape(IN_PATH % ('3_1_NES_ERM', )) 
smart_2 = readFilterShape(IN_PATH % ('3_1_ERM_NES', )) 

smart_3 = readFilterShape(IN_PATH % ('3_2_NES_OMM', )) 
smart_4 = readFilterShape(IN_PATH % ('3_2_OMM_NES', ))

smart_5 = readFilterShape(IN_PATH % ('3_3_ERM_OMM', )) 
smart_6 = readFilterShape(IN_PATH % ('3_3_OMM_ERM', )) 

smart_7 = readFilterShape(IN_PATH % ('3_4_NES_ERM_OMM', )) 
smart_8 = readFilterShape(IN_PATH % ('3_4_ERM_NES_OMM', ))
smart_9 = readFilterShape(IN_PATH % ('3_4_OMM_NES_ERM', ))


def get_site_shape_profile(SHAPE, Parser):
    import sys
    
    Start_Codon = []
    Stop_Codon = []
    for i in range(30):
        Start_Codon.append([])
        Stop_Codon.append([])
    for i in range(100):
        Start_Codon.append([])
        Stop_Codon.append([])
    
    for tid in SHAPE:
        try:
            feature = Parser.getTransFeature(tid)
        except:
            continue
        
        if feature['gene_type'] not in ('protein_coding', 'mRNA'):
            continue
        
        CDS_start = feature['cds_start']
        CDS_end = feature['cds_end']
        if CDS_start <= 50: continue                   
        if CDS_end - CDS_start <= 200: continue
        if feature['trans_len'] - CDS_end <= 50: continue
        
        for idx in range(CDS_start-30, CDS_start+100):
            shape_v = SHAPE[tid][idx]
            if shape_v != 'NULL':
                pos = idx - (CDS_start-30)
                Start_Codon[pos].append( float(shape_v) )
        
        for idx in range(CDS_end-100, CDS_end+30):
            shape_v = SHAPE[tid][idx]
            if shape_v != 'NULL':
                pos = idx - (CDS_end-100)
                Stop_Codon[pos].append( float(shape_v) )
    
    for idx in range(len(Start_Codon)):
        sys.stdout.writelines("%s " % (len(Start_Codon[idx]), ))
        #if len(Start_Codon[idx])>0: # Revised by Yuhan
        Start_Codon[idx] = sum(Start_Codon[idx])/len(Start_Codon[idx])
    print("")
    
    for idx in range(len(Stop_Codon)):
        sys.stdout.writelines("%s " % (len(Stop_Codon[idx]), ))
        Stop_Codon[idx] = sum(Stop_Codon[idx])/len(Stop_Codon[idx])
    print("")
    
    return Start_Codon, Stop_Codon

def FFT_transform(x,y):
    '''
    '''
    #L = len (y)                        # 信号长度
    #N = int(np.power(2,np.ceil(np.log2(L))))    # 下一个最近二次幂
    #FFT_y1 = np.abs(fft(y*2,N))/L*2      # N点FFT 变化,但处于信号长度
    #Fre = np.arange(int(N/2))*L/N 
    #FFT_y1 = FFT_y1[:int(N/2)]      # 取一半
    #return Fre, FFT_y1
    
    fft_y=fft(y)                          #快速傅里叶变换
    half_x = x[:int(len(x)//2)]  #取一半区间
    
    abs_y=np.abs(fft_y)                #取复数的绝对值，即复数的模(双边频谱)
    normalization_y=abs_y/len(x)            #归一化处理（双边频谱）                              
    normalization_half_y = normalization_y[:int(len(x)//2)]
    fwbest = normalization_half_y[signal.argrelextrema(normalization_half_y, np.greater)]
    xwbest = signal.argrelextrema(normalization_half_y, np.greater)
    return half_x,normalization_half_y*2



start_C_1, stop_C_1 = get_site_shape_profile(smart_1, hg38_parser) #ERM_LMNA
start_C_2, stop_C_2 = get_site_shape_profile(smart_2, hg38_parser) #LMNA_ERM

start_C_3, stop_C_3 = get_site_shape_profile(smart_3, hg38_parser) #ERM_OMM
start_C_4, stop_C_4 = get_site_shape_profile(smart_4, hg38_parser) #OMM_ERM

start_C_5, stop_C_5 = get_site_shape_profile(smart_5, hg38_parser) #LMNA_OMM
start_C_6, stop_C_6 = get_site_shape_profile(smart_6, hg38_parser) #OMM_LMNA

start_C_7, stop_C_7 = get_site_shape_profile(smart_7, hg38_parser) #ERM
start_C_8, stop_C_8 = get_site_shape_profile(smart_8, hg38_parser) #OMM
start_C_9, stop_C_9 = get_site_shape_profile(smart_9, hg38_parser) #LMNA

plt.figure(figsize=(15,10))

plt.subplot(4, 2, 1)
l1, = plt.plot(start_C_1, '-', color="#1f77b4")
l2, = plt.plot(start_C_2, '-', color="#ff7f0e")

plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
#plt.axvline(x=0, ymin=0, ymax = 1, linewidth=1, color='k')
#plt.axvline(x=1, ymin=0, ymax = 1, linewidth=1, color='k')
plt.axvline(x=29, ymin=0, ymax = 1, linewidth=2, color='k')
#plt.axvline(x=100, ymin=0, ymax = 1, linewidth=1, color='k')
plt.xticks([])
plt.title('NES vs ERM (-30bp)5UTR-CDS(+100bp) (4,786)')
plt.ylim(0.1, 0.5)



plt.subplot(4, 2, 2)
l3, = plt.plot(stop_C_1, '-', color="#1f77b4")
l4, = plt.plot(stop_C_2, '-', color="#ff7f0e")
plt.legend(handles = [l3, l4], labels = ['NES', 'ERM'], loc = 'best')

#plt.axvline(x=30, ymin=0, ymax = 1, linewidth=1, color='k')
#plt.axvline(x=129, ymin=0, ymax = 1, linewidth=1, color='k')
plt.axvline(x=99, ymin=0, ymax = 1, linewidth=2, color='k')
plt.xticks([])
plt.title('NES vs ERM (-100bp)CDS-3UTR(+30bp) (4,786)')
plt.ylim(0.1, 0.5)







plt.subplot(4, 2, 3)
l5, = plt.plot(start_C_3, '-', color="#1f77b4")
l6, = plt.plot(start_C_4, '-', color="#2ca02c")
plt.legend(handles = [l5, l6], labels = ['NES', 'OMM'], loc = 'best')

plt.axvline(x=29, ymin=0, ymax = 1, linewidth=2, color='k')
#plt.axvline(x=100, ymin=0, ymax = 1, linewidth=1, color='k')
plt.xticks([])
plt.title('NES vs OMM (-30bp)5UTR-CDS(+100bp) (4,216)')
plt.ylim(0.1, 0.6)



plt.subplot(4, 2, 4)
l7, = plt.plot(stop_C_3, '-', color="#1f77b4")
l8, = plt.plot(stop_C_4, '-', color="#2ca02c")
plt.legend(handles = [l7, l8], labels = ['NES', 'OMM'], loc = 'best')

#plt.axvline(x=29, ymin=0, ymax = 1, linewidth=1, color='k')
plt.axvline(x=99, ymin=0, ymax = 1, linewidth=2, color='k')
plt.xticks([])
plt.title('NES vs OMM (-100bp)CDS-3UTR(+30bp) (4,216)')
plt.ylim(0.1, 0.6)




plt.subplot(4, 2, 5)
l9, = plt.plot(start_C_5, '-', color="#ff7f0e")
l10, = plt.plot(start_C_6, '-', color="#2ca02c")
plt.legend(handles = [l9, l10], labels = ['ERM', 'OMM'], loc = 'best')

plt.axvline(x=29, ymin=0, ymax = 1, linewidth=2, color='k')
#plt.axvline(x=100, ymin=0, ymax = 1, linewidth=1, color='k')
plt.xticks([])
plt.title('ERM vs OMM (-30bp)5UTR-CDS(+100bp) (3,974)')
plt.ylim(0.1, 0.5)



plt.subplot(4, 2, 6)
l11, = plt.plot(stop_C_5, '-', color="#ff7f0e")
l12, = plt.plot(stop_C_6, '-', color="#2ca02c")
plt.legend(handles = [l11, l12], labels = ['ERM', 'OMM'], loc = 'best')

#plt.axvline(x=29, ymin=0, ymax = 1, linewidth=1, color='k')
plt.axvline(x=99, ymin=0, ymax = 1, linewidth=2, color='k')
plt.xticks([])
plt.title('ERM vs OMM (-100bp)CDS-3UTR(+30bp) (3,974)')
plt.ylim(0.1, 0.5)








plt.subplot(4, 2, 7)
l13, = plt.plot(start_C_7, '-', color="#1f77b4")
l14, = plt.plot(start_C_8, '-', color="#ff7f0e")
l15, = plt.plot(start_C_9, '-', color="#2ca02c")
plt.legend(handles = [l13, l14, l15], labels = ['NES', 'ERM', 'OMM'], loc = 'best')

plt.axvline(x=29, ymin=0, ymax = 1, linewidth=2, color='k')
#plt.axvline(x=100, ymin=0, ymax = 1, linewidth=1, color='k')
plt.xticks([])
plt.title('NES vs ERM vs OMM (-30bp)5UTR-CDS(+100bp) (3,206)')
plt.ylim(0.1, 0.5)



plt.subplot(4, 2, 8)
l16, = plt.plot(stop_C_7, '-', color="#1f77b4")
l17, = plt.plot(stop_C_8, '-', color="#ff7f0e")
l18, = plt.plot(stop_C_9, '-', color="#2ca02c")
plt.legend(handles = [l16, l17, l18], labels = ['NES', 'ERM', 'OMM'], loc = 'best')

#plt.axvline(x=29, ymin=0, ymax = 1, linewidth=1, color='k')
plt.axvline(x=99, ymin=0, ymax = 1, linewidth=2, color='k')
plt.xticks([])
plt.title('NES vs ERM vs OMM (-100bp)CDS-3UTR(+30bp) (3,206)')
plt.ylim(0.1, 0.5)

plt.tight_layout()
plt.savefig("./4_figure_origin.pdf")
plt.close()









plt.figure(figsize=(15,10))
shape_seq = [start_C_1,start_C_2,stop_C_1,stop_C_2,\
    start_C_3,start_C_4,stop_C_3,stop_C_4,\
    start_C_5,start_C_6,stop_C_5,stop_C_6,\
    start_C_7,start_C_8,start_C_9,stop_C_7,stop_C_8,stop_C_9]

plt.subplot(4, 4, 1)
x1,y1 = FFT_transform(list(range(30)),start_C_1[:30])
x2,y2 = FFT_transform(list(range(30)),start_C_2[:30])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs ERM (-30bp)5UTR (4,786)')

plt.subplot(4, 4, 2)
x1,y1 = FFT_transform(list(range(100)),start_C_1[30:])
x2,y2 = FFT_transform(list(range(100)),start_C_2[30:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs ERM CDS(+100bp) (4,786)')

plt.subplot(4, 4, 3)
x1,y1 = FFT_transform(list(range(100)),stop_C_1[:100])
x2,y2 = FFT_transform(list(range(100)),stop_C_2[:100])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs ERM (-100bp)CDS (4,786)')

plt.subplot(4, 4, 4)
x1,y1 = FFT_transform(list(range(30)),stop_C_1[100:])
x2,y2 = FFT_transform(list(range(30)),stop_C_2[100:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs ERM 3UTR(+30bp) (4,786)')

plt.subplot(4, 4, 5)
x1,y1 = FFT_transform(list(range(30)),start_C_3[:30])
x2,y2 = FFT_transform(list(range(30)),start_C_4[:30])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs OMM (-30bp)5UTR (4,216)')

plt.subplot(4, 4, 6)
x1,y1 = FFT_transform(list(range(100)),start_C_3[30:])
x2,y2 = FFT_transform(list(range(100)),start_C_4[30:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs OMM CDS(+100bp) (4,216)')

plt.subplot(4, 4, 7)
x1,y1 = FFT_transform(list(range(100)),stop_C_3[:100])
x2,y2 = FFT_transform(list(range(100)),stop_C_4[:100])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs OMM (-100bp)CDS (4,216)')

plt.subplot(4, 4, 8)
x1,y1 = FFT_transform(list(range(30)),stop_C_3[100:])
x2,y2 = FFT_transform(list(range(30)),stop_C_4[100:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('NES vs OMM 3UTR(+30bp) (4,216)')

plt.subplot(4, 4, 9)
x1,y1 = FFT_transform(list(range(30)),start_C_5[:30])
x2,y2 = FFT_transform(list(range(30)),start_C_6[:30])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('ERM vs OMM (-30bp)5UTR (3,974)')

plt.subplot(4, 4, 10)
x1,y1 = FFT_transform(list(range(100)),start_C_5[30:])
x2,y2 = FFT_transform(list(range(100)),start_C_6[30:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('ERM vs OMM CDS(+100bp) (3,974)')

plt.subplot(4, 4, 11)
x1,y1 = FFT_transform(list(range(100)),stop_C_5[:100])
x2,y2 = FFT_transform(list(range(100)),stop_C_6[:100])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('ERM vs OMM (-100bp)CDS (3,974)')

plt.subplot(4, 4, 12)
x1,y1 = FFT_transform(list(range(30)),stop_C_5[100:])
x2,y2 = FFT_transform(list(range(30)),stop_C_6[100:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
plt.legend(handles = [l1, l2], labels = ['NES', 'ERM'], loc = 'best')
plt.title('ERM vs OMM 3UTR(+30bp) (3,974)')

plt.subplot(4, 4, 13)
x1,y1 = FFT_transform(list(range(30)),start_C_7[:30])
x2,y2 = FFT_transform(list(range(30)),start_C_8[:30])
x3,y3 = FFT_transform(list(range(30)),start_C_9[:30])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
l3, = plt.plot(x3[2:],y3[2:], '-', color="#2ca02c")
plt.legend(handles = [l1, l2, l3], labels = ['NES', 'ERM', 'OMM'], loc = 'best')
plt.title('NES vs ERM vs OMM (-30bp)5UTR (3,206)')

plt.subplot(4, 4, 14)
x1,y1 = FFT_transform(list(range(100)),start_C_7[30:])
x2,y2 = FFT_transform(list(range(100)),start_C_8[30:])
x3,y3 = FFT_transform(list(range(100)),start_C_9[30:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
l3, = plt.plot(x3[2:],y3[2:], '-', color="#2ca02c")
plt.legend(handles = [l1, l2, l3], labels = ['NES', 'ERM', 'OMM'], loc = 'best')
plt.title('NES vs ERM vs OMM CDS(+100bp) (3,206)')

plt.subplot(4, 4, 15)
x1,y1 = FFT_transform(list(range(100)),stop_C_7[:100])
x2,y2 = FFT_transform(list(range(100)),stop_C_8[:100])
x3,y3 = FFT_transform(list(range(100)),stop_C_9[:100])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
l3, = plt.plot(x3[2:],y3[2:], '-', color="#2ca02c")
plt.legend(handles = [l1, l2, l3], labels = ['NES', 'ERM', 'OMM'], loc = 'best')
plt.title('NES vs ERM vs OMM (-100bp)CDS (3,206)')

plt.subplot(4, 4, 16)
x1,y1 = FFT_transform(list(range(30)),stop_C_7[100:])
x2,y2 = FFT_transform(list(range(30)),stop_C_8[100:])
x3,y3 = FFT_transform(list(range(30)),stop_C_9[100:])
l1, = plt.plot(x1[2:],y1[2:], '-', color="#1f77b4")
l2, = plt.plot(x2[2:],y2[2:], '-', color="#ff7f0e")
l3, = plt.plot(x3[2:],y3[2:], '-', color="#2ca02c")
plt.legend(handles = [l1, l2, l3], labels = ['NES', 'ERM', 'OMM'], loc = 'best')
plt.title('NES vs ERM vs OMM 5UTR(+30bp) (3,206)')



plt.tight_layout()
plt.savefig("./4_figure_FFT.pdf")
plt.close()


