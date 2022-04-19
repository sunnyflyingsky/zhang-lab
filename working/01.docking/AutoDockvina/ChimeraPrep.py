from chimera import runCommand as rc
from MolKit import Read
import MolKit.molecule
import MolKit.protein
import chimera
import DockPrep
import AddH
from AddCharge import estimateFormalCharge, addNonstandardResCharges
from DockPrep import prep
from WriteDMS import writeDMS

import os
import sys
import getopt

def prep_receptor(name,inputPath,outputPath):
    '''
    '''
    rc("open "+inputPath+'/'+name+".pdb") #
    rc("delete ligand")
    rc("delete ions")
    mols = chimera.openModels.list(modelTypes=[chimera.Molecule])
    if len(mols)>1:
        rc('delete ~#0.1')
    else:
        pass
    receptor = chimera.openModels.list(modelTypes=[chimera.Molecule])
    prep(receptor, method='gas', addHFunc=AddH.simpleAddHydrogens)
    rec_charge = estimateFormalCharge(receptor[0].atoms)
    l,w,h,center = calculateWindows(receptor[0].atoms)
    rc("write format pdb #0 "+outputPath+'/'+name+"_receptor.pdb")
    rc("delete H")
    rc("surface #0")
    surf = chimera.openModels.list(modelTypes=[chimera.MSMSModel])[0]
    writeDMS(surf, outputPath+"/receptor.ms")
    rc("close all")
    return l,w,h,center


def prep_ligand(name,inputPath,outputPath):
    '''
    '''
    rc("open "+inputPath+'/'+name+".pdb") #
    rc("select ligand")
    rc("delete ~selected")
    mols = chimera.openModels.list(modelTypes=[chimera.Molecule])
    if len(mols)>1:
        rc('delete ~#0.1')
    else:
        pass
    ligand = chimera.openModels.list(modelTypes=[chimera.Molecule])
    residues = ligand[0].residues
    if len(residues)>1:
        rc('delete ~#0:.A')
        ligand = chimera.openModels.list(modelTypes=[chimera.Molecule])
    else:
        pass
    prep(ligand, method='gas', addHFunc=AddH.simpleAddHydrogens)
    lig_charge = estimateFormalCharge(ligand[0].atoms)
    rc("write format pdb #0 "+outputPath+'/'+name+"_ligand.pdb")
    rc("close all")

def calculateWindows(atomLists):
    a0 = atomLists[0]
    x,y,z = a0.coord()[0],a0.coord()[1],a0.coord()[2]
    scale_x = [x,x]
    scale_y = [y,y]
    scale_z = [z,z]
    center = [0,0,0]
    for a in atomLists:
        x,y,z = a.coord()[0],a.coord()[1],a.coord()[2]
        if x < scale_x[0]:
            scale_x[0]=x
        elif x > scale_x[1]:
            scale_x[1]=x
        if y < scale_y[0]:
            scale_y[0]=y
        elif y > scale_y[1]:
            scale_y[1]=y
        if z < scale_z[0]:
            scale_z[0]=z
        elif z > scale_z[1]:
            scale_z[1]=z
        center[0]+=x
        center[1]+=y
        center[2]+=z
    for i in range(len(center)):
        center[i] /= len(atomLists)
    
    return (scale_x[1]-scale_x[0])*1.2+4.0,(scale_y[1]-scale_y[0])*1.2+4.0, \
        (scale_z[1]-scale_z[0])*1.2+4.0,center


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
    name = info[-1][:-4]
    prep_ligand(name,inputPath,outputPath)
    print('get ligand!')
    l,w,h,center = prep_receptor(name,inputPath,outputPath)
    with open(outputPath+'/windows.txt','w') as write_object:
        write_object.write('x_lenghth: {}\n'.format(l))
        write_object.write('y_width: {}\n'.format(w))
        write_object.write('z_height: {}\n'.format(h))
        write_object.write('center: {} {} {}\n'.format(center[0],center[1],center[2]))
    print('get receptor!')
