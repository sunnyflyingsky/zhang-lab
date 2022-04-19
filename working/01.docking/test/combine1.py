import os
import chimera
import DockPrep
import sys
from DockPrep import prep
from WriteDMS import writeDMS
from WriteMol2 import writeMol2
import Midas
from chimera import runCommand as rc
import AddH
#chimera --nogui --nostatus --script "combine1.py 1_Rec_noH.pdb 3_predictions_scored_best.pdb"

Rec = "1_Rec_noH.pdb"
Lig = "3_predictions_scored_best.pdb"

rc("open "+str(Rec))
rc("open "+str(Lig)) 

Rec = chimera.openModels.list(modelTypes=[chimera.Molecule])
Midas.write(Rec,None,"4_results_tmp.pdb")

