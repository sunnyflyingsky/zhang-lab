from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw

drugs = {}
with open('/data/zhangjiasheng/data/shape2structure/04_RNA_drug/drug2smile.txt') as read_object:
    for line in read_object:
        info = line.strip().split('\t')
        if len(info)==3:
            try:
                mol = Chem.MolFromSmiles(info[1])
            except:
                mol = Chem.MolFromSmiles(info[2])

        elif len(info)==2:
            mol = Chem.MolFromSmiles(info[1])
        if mol == None:
            print(info[0])
            continue
        smile = Chem.MolToSmiles(mol)
        drugs[info[0]]=smile
        #print(info[1]==smile)
with open('/data/zhangjiasheng/toolkits/RNALigands/data/miRBase.txt') as read_object:
    for line in read_object:
        info = line.strip().split('\t')
        name = info[-1]
        mol = Chem.MolFromSmiles(name)
        smile = Chem.MolToSmiles(mol)
        if name not in drugs.keys():
            drugs[name]=smile
with open('/data/zhangjiasheng/data/shape2structure/04_RNA_drug/drug2smile_corrected.txt','w') as write_object:
    for key,value in drugs.items():
        write_object.write('{}\t{}\n'.format(key,value))








