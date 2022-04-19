#!/usr/bin/python
#python 4_fix.py -i 3_SS_revised_details.txt > 4_SS_revised_details_removed.txt
import sys
import re
for arg in range(len(sys.argv)):
    if sys.argv[arg]=='-i':
        data=sys.argv[arg+1]



def read_structure(filename):
    ids, strands, seqs, structures= [], [], [], []
    for line in open(filename,'r'):
            id, strand, seq, structure= line.split('\t')
            ids.append(id)
            strands.append(strand)
            seqs.append(seq)
            structures.append(structure.replace('\n',''))
    return ids, strands, seqs, structures



s_ids, s_strands, s_seqs, s_structures= read_structure(data)



#print chromosome
#w=open(Out,'w')
for i in range(0,len(open(data).readlines()),1):
    #if str(s_seqs[i][-14:])=="&.&.&.&.&.&.&.":
    if (s_ids[i]=="2kd4" or s_ids[i]=="5eme" or s_ids[i]=="5emf" or s_ids[i]=="5nz6"):
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i].replace('&','') + '\t' + s_structures[i].replace('&',''))
    elif (s_ids[i]=="1e8s"):
        seq="UCgcuagcgagaccccgucucugccgggcgcgguggcgcgcgccuguagucccagcuacucgggaggcugaggugggaggaucgcuag"
        dbn="--......................................................................................"
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + seq + '\t' + dbn)
    elif (s_ids[i]=="5lmo"):
         print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-5] + '\t' + s_structures[i][:-5])
    elif re.match("&.&.&.&.", s_seqs[i][-8:]):
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-8] + '\t' + s_structures[i][:-8])
    elif re.match("&.&.&.", s_seqs[i][-6:]):
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-6] + '\t' + s_structures[i][:-6])
    elif re.match("&.&.", s_seqs[i][-4:]):
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-4] + '\t' + s_structures[i][:-4])
    elif re.match("&.", s_seqs[i][-2:]):
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-2] + '\t' + s_structures[i][:-2])
    #elif str(s_seqs[i][-3:])=="&AA": #5LMO
    #    print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-3] + '\t' + s_structures[i][:-3])        
    #elif str(s_seqs[i][-3:])=="Cc&": #5NZ6
    #    print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-3] + '\t' + s_structures[i][:-3])       
    else:
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i] + '\t' + s_structures[i])

    """
    elif str(s_seqs[i][-12:])=="&.&.&.&.&.&.":
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-12] + '\t' + s_structures[i][:-12])
    elif str(s_seqs[i][-10:])=="&.&.&.&.&.":
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-10] + '\t' + s_structures[i][:-10])
    elif str(s_seqs[i][-8:])=="&.&.&.&.":
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-8] + '\t' + s_structures[i][:-8])
    elif str(s_seqs[i][-6:])=="&.&.&.":
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-6] + '\t' + s_structures[i][:-6])
    elif str(s_seqs[i][-4:])=="&.&.&.":
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-4] + '\t' + s_structures[i][:-4])
    elif str(s_seqs[i][-2:])=="&.&.":
        print(s_ids[i]+ '\t' + s_strands[i] + '\t' + s_seqs[i][:-2] + '\t' + s_structures[i][:-2])

        
    #else:
    #    print("error_unknown")
    """