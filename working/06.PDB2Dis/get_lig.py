import os
import sys
import atomium
#python Lig_prep2.py 3f2t.pdb 11

metals = ["LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE",
          "CO", "NI", "CU", "ZN", "HA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU",
          "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA", "LA", "CE", "PR", "ND", "PM",
          "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W",
          "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "FR", "RA", "AC",
          "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO",
          "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "UUT", "FL", "LV"]

others = ["WO2", "OHX"]

data = sys.argv[1]
cutoff = sys.argv[2]

pdb = atomium.open(str(data))

ligands=list(pdb.model.ligands())

ligands_new=[]
strands_new=[]
for i in range(len(ligands)):
    ligand=str(ligands[i]).split(" ")[0]
    strand=str(ligands[i]).split(" ")[-1].split(".")[0]
    new=ligand+"_"+strand
    ligands_new.append(new)
    strands_new.append(strand)

ligands_clean=[]
for i in range(len(ligands_new)):
    #sm_name=str(list(pdb.model.ligands())[i]).split(" ")[0]
    sm=ligands_new[i]
    sm_name=str(sm).split("_")[0]
    #print(sm_name)
    if sm_name not in metals and sm_name not in others:
        ligands_clean.append(sm)
        
strands_new = list(dict.fromkeys(strands_new))
ligands_real=list(filter(lambda s: ligands_clean.count(s) <= int(cutoff), ligands_clean))
#print(ligands_real)
#print(strands_new)
Ligs=list(set(ligands_real))
print(Ligs)