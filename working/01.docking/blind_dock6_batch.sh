#source /Share2/home/zhangqf7/yuhan/env.sh && conda activate dock
OUT=/Share2/home/zhangqf7/yuhan/6_Docking/dock6/results_pdb
TOOLS=${OUT}/tools
INI=${OUT}/ini

DATA=/Share2/home/zhangqf7/yuhan/6_Docking/20220225_datasets_PDB/$1
DRUG=/Share2/home/zhangqf7/yuhan/6_Docking/20220225_datasets_PDB/Drug_pdb

file=$2
filename=$(echo $file | cut -d . -f1)
final_data=${OUT}/$3

cd ${OUT}
mkdir ${final_data}
chimera --nogui --script "${TOOLS}/dockprep.py ${DATA}/${file} ${OUT}"
cd ${OUT}/$filename
sh ${TOOLS}/grid.sh
for ligand in `ls ${DRUG} --color=never | sort -V`
do
ligname=$(echo $ligand | cut -d . -f1)

mkdir ${OUT}/$filename/Ligands/${ligand}
cd ${OUT}/$filename/Ligands/${ligand}

chimera --nogui --script "${TOOLS}/ligprep.py ${DRUG}/${ligand} ${OUT}"

#timeout 120m bash ${TOOLS}/docking.sh ${INI}/flexiable.in
timeout 180m bash ${TOOLS}/docking.sh ${INI}/flexiable.in


chimera --nogui --script "${TOOLS}/rmH.py 1_Lig.pdb 1_Lig_noH.pdb"
chimera --nogui --script "${TOOLS}/rmH.py 3_predictions_scored.pdb 3_predictions_scored_noH.pdb"

awk -v na=3_predictions_scored 'BEGIN{id=1}/^MODEL/{id=$2;}{print $0 > na"_"id".pdb" ;}' 3_predictions_scored_noH".pdb"

cp 3_predictions_scored_1.pdb 3_predictions_scored_best.pdb

weight=$(obprop 3_predictions_scored_best.pdb | grep ^mol_weight | awk '{print $2}')
#encode=$(obprop 3_predictions_scored_best.pdb | grep ^canonical_SMILES | awk '{print $2}')

#binding affinity
grep Grid_Score 3_predictions_scored.mol2 | awk '{print $3}' | head -n 1 | awk -v rec_name="${filename}" -v lig_name="${ligand}" -v lig_weight="${weight}" '{print rec_name "\t" lig_name "\t" lig_weight "\t" $1 }' > 3_predictions_score_best.txt

chimera --nogui --nostatus --script "${TOOLS}/combine.py ${OUT}/$filename/1_Rec_noH.pdb 3_predictions_scored_best.pdb ${TOOLS}" #4_results_final.pdb
chimera --nogui --nostatus --script "${TOOLS}/rmH.py 4_results_final.pdb 4_results_final_noH.pdb"
chimera --nogui --nostatus --script "${TOOLS}/combine2.py ${OUT}/$filename/1_Rec_noH.pdb 3_predictions_scored_best.pdb" #4_results_tmp.pdb

python ${TOOLS}/merge.py -i 4_results_final_noH.pdb -r 4_results_tmp.pdb -o 5_results.pdb

score=$(grep Grid_Score 3_predictions_scored.mol2 | awk '{print $3}' | head -n 1 | bc -l | xargs printf "%.1f")

cp 5_results.pdb ${final_data}/${filename}_${ligname}_${score}.pdb

cat ${OUT}/$filename/Ligands/${ligand}/3_predictions_score_best.txt >> ${OUT}/$filename/3_predictions_score_best_all.txt
done

