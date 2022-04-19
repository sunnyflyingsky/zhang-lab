# python sdftosmiles.py molecules.sdf

import sys
from rdkit import Chem

import warnings
warnings.filterwarnings("error")

def converter(file_name):
    mols = [ mol for mol in Chem.SDMolSupplier( file_name ) ]
    #outname = file_name.split(".sdf")[0] + ".smi"
    #out_file = open( outname, "w" )
    for mol in mols:
        try:
            smi = Chem.MolToSmiles(mol,isomericSmiles=False) #https://www.jianshu.com/p/c0df2942d8d1
            name = mol.GetProp("_Name")
        except:
            print( "{}\t{}\t{}".format(file_name.split("_")[0].split("/")[-1],"null","error"))
        else:
            #out_file.write( "{}\t{}\n".format(smi, name ))
            print( "{}\t{}\t{}".format(name,smi,"OK"))
    #out_file.close()

if __name__=="__main__":
    converter( sys.argv[1] )