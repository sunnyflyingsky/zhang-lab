Dock6, autodock vina, MORDOR: 
PDB_id: 1aju, 2esj, 3f2t


/Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking
/home/yuhan/zhangjiasheng/dynamic_region_analysis_NES_ERM

如果想加环境变量的话，就加在~/yuhan/env.sh
source ~/yuhan/env.sh
conda activate dock

scp -r /Share2/home/zhangqf7/yuhan/0_genome/STAR/mm10_Gencode yuhan@20.20.91.14:/home/yuhan/zhaokang_SNP/ref
rsync -av --progress /path/to/dir dongzhuoer@20.20.91.18:/data1/raw/dir （默认密码与用户名相同）

v2 start
==================================

bqueues
bjobs
bhosts
bkill

501-505

522-540

bsub -o output.txt -e err.log -q Z-ZQF -m node533 -n 1 "xxx"

bsub -q Z-ZQF -m node533 -n 1 "xxx"


bsub -q Z-BNODE -n 1 "xxx"
bsub -q Z-HNODE -n 1 "xxx"

===================================

https://github.com/sahrendt0/Scripts/blob/master/docking_scripts/prepare_receptor4.py
pythonsh xxx/prepare_receptor4.py -r xxx.pdb -o xxx.pdbqt

chimera --nogui --script "xxx.py"
在chimera的IDLE里面，可以通过from chimera import runCommand as rc,去执行一些操作，比如rc("select xxx")

babel -i pdbqt xxx.pdbqt -o pdb xxx.pdb

import chimera
import DockPrep
from DockPrep import prep
Ligand = chimera.openModels.list(modelTypes=[chimera.Molecule])
prep(Ligand, method='gas', addHFunc=AddH.simpleAddHydrogens)

软件计算化学ChemAxon
能量量化计算，实验参考。datamining。

=====================================
GPU
login: ssh -Y -p 2014 yuhan@101.6.120.41 -L 127.0.0.1:1111:127.0.0.1:1111
password:feiyuhan
Juypter: jupyter-notebook --no-browser --port 1111

显示显存占用
nvidia-smi dmon -s p -i 0 -c 20 |grep -v "^#"|awk 'BEGIN{s=0}{s+=$2}END{print( s/NR "W")}'  
nvidia-smi dmon -s u -i 0 -c 20 |grep -v "^#"|awk 'BEGIN{s=0}{s+=$2}END{print( s/NR "%")}'






