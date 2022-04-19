#!/usr/bin/python
import sys
#python 3_merge.py -i SS_revised.txt -j 2_SS_details.txt > SS_revised_details.txt

for arg in range(len(sys.argv)):
    if sys.argv[arg]=='-i':
        data=sys.argv[arg+1]
    elif sys.argv[arg]=='-j':
        detail=sys.argv[arg+1]


def read_files(filename):
    ids, strands, seqs, structures = [], [], [], []
    for line in open(filename,'r'):
        id, strand, seq, structure = line.split('\t')
        ids.append(id)
        strands.append(strand)
        seqs.append(seq)
        structures.append(structure.replace('\n',''))

    return ids, strands, seqs, structures


def read_details(filename):
    ids_d, strands_d, seqs_d, structures_d= [], [], [], []
    for line in open(filename,'r'):
            id, strand, seq, structure= line.split('\t')
            ids_d.append(id)
            strands_d.append(strand)
            seqs_d.append(seq)
            structures_d.append(structure.replace('\n',''))
    return ids_d, strands_d, seqs_d, structures_d



def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

ids, strands, seqs, structures = read_files(data)


ids_d, strands_d, seqs_d, structures_d = read_details(detail)



#print chromosome
#w=open(Out,'w')
for i in range(0,len(open(data).readlines()),1):
    for j in range(0,len(open(detail).readlines()),1):
        if str(ids[i])==str(ids_d[j]) and str(strands[i]) == str(strands_d[j]) :
            #print(ids[i]+ '\t' + nums[i] + '\t' + cats[i] + '\t' + species[i] + '\t' + resolutions[i] + '\t' + dates[i] + '\t' + ligands[i])
            l=list(find_all(seqs[i],'&'))
            if l:
                for k in l:
                    seqs_d[j] = seqs_d[j][:k] + "&" + seqs_d[j][k:]
                    structures_d[j] = structures_d[j][:k] + "&" + structures_d[j][k:] 
                print(ids[i]+"\t"+strands[i]+"\t"+seqs[i]+"\t"+structures_d[j])
            else:
                print(ids[i]+"\t"+strands[i]+"\t"+seqs[i]+"\t"+structures_d[j])

        #else:
        #    print("error_unknown")