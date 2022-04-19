from chimera import runCommand as rc
from MolKit import Read
import MolKit.molecule
import MolKit.protein
import chimera
import DockPrep
import AddH
from AddCharge import estimateFormalCharge, addNonstandardResCharges
from DockPrep import prep
from chimera import match

import os
import sys
import getopt


def ligand_native(name,inputPath,outputPath):
    '''
    '''
    print(inputPath+"/"+name+"_ligand.pdb")
    rc("open "+inputPath+"/"+name+"_ligand.pdb")
    rc("delete H")
    rc("write format pdb #0 "+outputPath+'/'+name+"_ligand.pdb")
    rc("close all")

def ligand_dock(name,inputPath,outputPath):
    '''
    '''
    rc("open "+outputPath+'/'+name+".pdb")
    #rc("delete H")
    mols = chimera.openModels.list(modelTypes=[chimera.Molecule])
    for i,mol in enumerate(mols):
        num = str(i+1)
        rc("write format pdb #0."+num+" "+outputPath+'/'+name+"_ligand_d"+num+".pdb")
    rc("close all")
    return int(num)

def cal_rmsd(name,inputPath,outputPath,num):
    write_object = open(outputPath+"/RMSD_summary.txt",'w')
    for i in range(num):
        rc("open "+inputPath+"/"+name+"_ligand.pdb")
        num = str(i+1)
        rc("open "+inputPath+"/"+name+"_ligand_d"+str(num)+".pdb")
        rc("delete H")
        #rc("rmsd #0 #1")
        ligand = chimera.openModels.list(modelTypes=[chimera.Molecule])
        a = ligand[0].atoms
        b = ligand[1].atoms
        t = match.matchAtoms(a,b)
        write_object.write('d{}: {}\n'.format(num,t[1]))
        rc('close all')
    write_object.close()


Usage = 'Usage: ' + sys.argv[0]
Usage += '''
    <Requires>
    -i input pdb file path
    -o output path

    [Options]
    None
'''
Usage += "EX: " + sys.argv[0] + ''' -i /Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/data_set/xxxx.pdb 
    /Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/temp_set'''
if len(sys.argv)<5 or not sys.argv[1].startswith('-'):sys.exit(Usage)

if __name__ =='__main__':
    print('run!')
    oplist,alist = getopt.getopt(sys.argv[1:],'hi:o:')
    for opt in oplist:
        if opt[0] == '-h':sys.exit(Usage)
        elif opt[0] == '-i':inputInfo = opt[1]
        elif opt[0] == '-o':outputPath = opt[1]
        else:
            sys.exit(Usage)
    info = inputInfo.strip().split('/')
    inputPath = '/'.join(info[:-1])
    #print(inputPath)
    name = info[-1][:-4]
    #ligand_native(name,inputPath,outputPath)
    num = ligand_dock(name,inputPath,outputPath)
    cal_rmsd(name,outputPath,outputPath,num)
