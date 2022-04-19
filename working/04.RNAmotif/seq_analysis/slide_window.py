import sys
window=30
step=10
#python slide_window.py -i results_final_anno_3utr_merge.txt > results_final_anno_3utr_merge_window.txt
for arg in range(len(sys.argv)):
    if sys.argv[arg]=='-i':
        data=sys.argv[arg+1]
    if sys.argv[arg]=='-w':
        window=int(sys.argv[arg+1])
    if sys.argv[arg]=='-s':
        step=int(sys.argv[arg+1])

def read_files(filename):
    ids, starts, ends, dyn_nuls, val_nuls, avgs, labels, p_values = [], [], [], [], [], [], [], []
    read_object = open(filename,'r')
    read_object.readline()
    for line in read_object:
        id, start, end, dyn_nul, val_nul, avg, label, p_value = line.split('\t')
        if int(end)-int(start) >=20 and  int(end)-int(start)<30:
            label = label + '^'*(30-int(end)+int(start))
            end = str(int(start)+30)
        ids.append(id)
        starts.append(start)
        ends.append(end)
        dyn_nuls.append(dyn_nul)
        val_nuls.append(val_nul)
        avgs.append(avg)
        labels.append(label)
        p_values.append(p_value.replace('\n',''))
    return ids, starts, ends, dyn_nuls, val_nuls, avgs, labels, p_values


ids, starts, ends, dyn_nuls, val_nuls, avgs, labels, p_values = read_files(data)

print('id'+'\t'+'start'+'\t'+'end'+'\t'+'dyn_nul'+'\t'+'val_nul'+'\t'+'avg'+'\t'+'label'+'\t'+'p_value')
for i in range(0,len(ids),1):
    #
    for j in range(int(starts[i]),int(ends[i]) - int(window)+1,int(step)):
        #print(i,j)
        print(str(ids[i])+"\t"+str(j)+"\t"+str(int(j)+int(window))+"\t"+str(dyn_nuls[i])+\
            "\t"+str(val_nuls[i])+"\t"+str(avgs[i])+"\t"+str(labels[i][j-int(starts[i]):j-int(starts[i])+window])+"\t"+str(p_values[i]))