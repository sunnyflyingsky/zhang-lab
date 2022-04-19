import GAP
import General
hg38_parser = GAP.init("/home/yuhan/smartSHAPE_data/get_seq/hg38.genomeCoor.bed", "/home/yuhan/smartSHAPE_data/get_seq/hg38_transcriptome.fa")

window=101
step=70
file = open('/home/yuhan/smartSHAPE_data/Human_293T_smartSHAPE.out')
f = open("/home/yuhan/smartSHAPE_data/Human_293T_smartSHAPE_seq.txt", 'w')
for line in file:
    #print(line, end="")
    trans_id=line.strip()
    seq = hg38_parser.getTransSeq(trans_id)
    shape = General.load_shape("/home/yuhan/smartSHAPE_data/Human_293T_smartSHAPE.out")[trans_id]
    for i in range(0,len(seq) - window,step):
        if(len(seq)==len(shape)):
            a = seq[i:(i+window)]
            b = shape[i:(i+window)]
            if "NULL" not in b:
                b=','.join(b)
                #print(str(a)+"\t"+str(b))
                f.write(str(trans_id)+"\t"+str(a)+"\t"+str(b)+"\n")
        else:
            print("wrong_length")