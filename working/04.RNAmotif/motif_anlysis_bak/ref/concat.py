import os
import re

def concat():
    '''
    结合多个MEME文件成一整个文件
    ''' 
    refpath = 'RNAmotif/motif_anlysis/ref/oRNAment_human/RBP_id_encoding.csv'
    outpath1 = 'RNAmotif/motif_anlysis/ref/oRNAment_human/RBP_human.txt'
    outpath2 = 'RNAmotif/motif_anlysis/ref/oRNAment_human/Homo_sapiens.meme'
    inputpath = 'RNAmotif/motif_anlysis/ref/oRNAment_human/MEMEs'
    nameLists = os.listdir(inputpath)
    print(len(nameLists))
    write_object = open(outpath2,'w')
    flag = 0
    for name in nameLists:
        with open(inputpath+'/'+name) as read_object:
            if flag:
                for i in range(9):
                    read_object.readline()
            for line in read_object:
                if line.startswith('MOTIF'):
                    info = line.strip().split(' ')
                    info[1] = 'M'+str(name.replace('.meme',''))+'_0.8'
                    info[2] = '('+info[2]+')'
                    write_object.write('{}\n'.format(' '.join(info)))
                else:
                    write_object.write(line)
            flag = 1
    write_object.close()
    write_object = open(outpath1,'w')
    with open(refpath) as read_object:
        for line in read_object:
            info = line.strip().split(',')
            num = int(info[0])
            if num<=9:
                info[0] = 'M00'+str(num)+'_0.8'
            elif num<=99:
                info[0] = 'M0'+str(num)+'_0.8'
            else:
                info[0] = 'M'+str(num)+'_0.8'
            info[1] = info[1].strip().replace('"','').replace(' ','_')
            write_object.write('T0000\t{}\tCollapsed\tHuman\n'.format('\t'.join(info)))
    write_object.close()

if __name__ == '__main__':
    concat()
