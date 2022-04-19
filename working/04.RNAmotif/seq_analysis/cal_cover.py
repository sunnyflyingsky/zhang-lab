
def calc_cover(inpath,refpath,outpath):
    '''
    '''
    anno = {}
    with open(refpath) as read_object:
        for line in read_object:
            if line.startswith('@SQ'):
                info = line.strip().split('\t')
                anno[info[1].split(':')[1]] = int(info[2].split(':')[1])
    inlength = {}
    with open(inpath) as read_object:
        for line in read_object:
            info = line.strip().split('\t')
            l = int(info[2]) - int(info[1])
            try:
                inlength[info[0]]+=l
            except:
                inlength[info[0]]=l
    al,av = 0,0
    with open(outpath+'/coverage_all.txt','w') as write_object:
        for key,value in anno.items():
            if key in inlength.keys():
                c=round(inlength[key]/value,6)
                al+=inlength[key]
            else:
                c=0
                al+=0
            av+=value
            write_object.write('{}\t{}\n'.format(key,c))
        c = round(al/av,6)
        write_object.write('{}\t{}\n'.format('all',c))





if __name__ == '__main__':
    inpath = 'RNAmotif/seq_analysis/data/3_merged2_large_unique_sorted_merged.bed'
    refpath = 'RNAmotif/seq_analysis/data/3_merged2.fastq.sorted_bc1_anno.txt'
    outpath = 'RNAmotif/seq_analysis/res'
    calc_cover(inpath, refpath, outpath)