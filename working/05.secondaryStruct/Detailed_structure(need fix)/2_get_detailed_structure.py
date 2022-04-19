import sys

#cd /home/yuhan/dssr/2.3_SS_details
#for file in `ls | sort`; do echo ${file} ; python /home/yuhan/dssr/2_get_detailed_structure.py ${file} >> /home/yuhan/dssr/2_SS_details.txt; done


name=sys.argv[1]
#next(file)
#seq=[]
#ss=[]
strand_tmp=[]

file = open(name)
for line in file:
    line_re=line.strip("\n")
    line_re2="\t".join(line_re.split())
    element=line_re2.split("\t")
    strand_tmp.append(element[3].split(".")[0])
strand = list(set(strand_tmp))    
#print(strand)

#file2 = open(name)
for i in range(len(strand)):
    file2 = open(name)
    seq=[]
    ss=[]
    for line in file2:
        line_re=line.strip("\n")
        line_re2="\t".join(line_re.split())
        element=line_re2.split("\t")
        if element[3].split(".")[0]==strand[i]:
            if element[2]=="(":
                seq.append(element[1])
                ss.append(element[2])
            elif element[2]==")":
                seq.append(element[1])
                ss.append(element[2])
            elif(element[2]=="."):
                #print(element[0]+"\t"+element[1]+"\t"+element[5])
                if 'junction-loop' in element[5]:
                    #print(element[0]+"\t"+element[1]+"\t"+"J"+"\t"+"junction")  
                    seq.append(element[1])
                    ss.append("J")
                elif 'internal-loop' in element[5]:
                    #print(element[0]+"\t"+element[1]+"\t"+"I"+"\t"+"internal")
                    seq.append(element[1])
                    ss.append("I")
                elif 'hairpin-loop' in element[5]:
                    #print(element[0]+"\t"+element[1]+"\t"+"H"+"\t"+"hairpin")
                    seq.append(element[1])
                    ss.append("H")
                elif 'bulge' in element[5]:
                    #print(element[0]+"\t"+element[1]+"\t"+"L"+"\t"+"bulge")
                    seq.append(element[1])
                    ss.append("L")
                elif 'ss-non-loop' in element[5]:
                    #print(element[0]+"\t"+element[1]+"\t"+"."+"\t"+"unpair")  
                    seq.append(element[1])
                    ss.append(".")
    #print(name.split(".")[0],strand[i],"".join(seq),"".join(ss))
    print(name.split(".")[0]+"\t"+strand[i]+"\t"+"".join(seq)+"\t"+"".join(ss))

