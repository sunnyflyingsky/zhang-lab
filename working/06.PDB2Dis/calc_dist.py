from getSequence import cal,ext
from pdbReader import PDB
from settings import PDB_Entry,RingData
import math

def neiborSeq(inputpath:str='3DIR.txt'):
    '''
    '''
    file = open(inputpath)
    extend=7
    win=3
    cutoff=10
    remove=8
    #next(file)
    print("PDB_ID"+"\t"+"Ligand_name"+"\t"+"SMILES"+"\t"+"start"+"\t"+"end"+"\t"+"Sequence"+"\t"+"Structure"+"\t"+"RMSD"+"\t"+"Affinity"+"\t"+"Weight")
    for line in file:
        line_re=line.strip("\n")
        element=line_re.split("\t")
        #print(element[0]+"\t"+element[1]+"\t"+element[2]+"\t"+element[3]+"\t"+element[4]+"\t"+element[5]+"\t"+element[6])
        
        if (len(element[3]))>cutoff:
            cal(element,extend,win)
        elif (len(element[3]))<=cutoff and (len(element[3]))>=remove:
            ext(element,extend)
        elif (len(element[3]))<remove:
            pass

def calDist(AtomList:list,drugCenter:list):
    '''
    '''
    distance = []
    #for i in range(3):
    #    distance += (AtomList[0][i]-drugCenter[i])**2
    #distance = math.sqrt(distance)
    for atom in AtomList:
        d = 0
        for i in range(3):
            d += (atom[i]-drugCenter[i])**2
        d = math.sqrt(d)
        distance.append(d)
    #    if d<distance:
    #        distance = d
    distance = sorted(distance)
    #print(distance[:10])
    #print(distance[660:])
    return distance[len(AtomList)//2]



if __name__ == '__main__':
    #print('run!')
    neiborSeq('PDB2Dis/res/3DIR.txt')
    write_obejct = open('PDB2Dis/res/3DIR_seq_final.txt','w')
    with open('PDB2Dis/res/3DIR_seq.txt',encoding='utf-8') as read_object:
        line = read_object.readline()
        write_obejct.write(line.strip()+'\t'+'distance\n')
        for line in read_object:
            info = line.strip().split('\t')
            if len(info)<2:
                continue
            PDB_ID = info[0].strip()
            Ligand_name = info[1].strip()
            start = int(info[3].strip())
            end = int(info[4].strip())
            sequence = info[5].strip()
            drugCenter = [0,0,0]
            AtomList = []
            RNA = PDB('PDB2Dis/'+PDB_ID+'.pdb')
            Ligand = PDB('PDB2Dis/'+PDB_ID+'.pdb','ligand')
            for j in range(len(Ligand.residueList)):
                for key,value in Ligand.residueList[j].items():
                    if value==Ligand_name:
                        i = 0
                        for Ai in Ligand.Conformers_Atom[j][key].values():
                            drugCenter[0] += Ai.X
                            drugCenter[1] += Ai.Y
                            drugCenter[2] += Ai.Z
                            i+=1
                        drugCenter[0]/=i
                        drugCenter[1]/=i
                        drugCenter[2]/=i
                        chainId = j
                        break
            if chainId >= len(RNA.residueList):
                chainId = 0
            s0 = min(RNA.residueList[chainId].keys())
            for i in range(start+s0,end+s0,1):
                if i not in RNA.Conformers_Atom[chainId].keys():
                    continue
                for Ai in RNA.Conformers_Atom[chainId][i].values():
                    AtomList.append([0,0,0])
                    AtomList[-1][0] = Ai.X
                    AtomList[-1][1] = Ai.Y
                    AtomList[-1][2] = Ai.Z
            distance = calDist(AtomList,drugCenter)
            write_obejct.write(line.strip()+'\t'+str(distance)+'\n')


