cd /Share2/home/zhangqf7/yuhan/6_Docking/dock6/results_pdb
source /Share2/home/zhangqf7/yuhan/env.sh && conda activate dock




SOURCE=RNA_80
FINAL=final_80
DATA=/Share2/home/zhangqf7/yuhan/6_Docking/20220225_datasets_PDB/${SOURCE}
for file in `ls ${DATA} --color=never | sort -n`
do
echo ${file}
bsub -q Z-ZQF -m node501 -n 1 "sh blind_dock6_batch.sh ${SOURCE} ${file} ${FINAL}"
done