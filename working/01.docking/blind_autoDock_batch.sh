OUT=/Share2/home/zhangqf7/yuhan/6_Docking/autodock/results_pdb
TOOLS=${OUT}/tools
INI=${OUT}/ini


Utilities=${OUT}/programs/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24


DATA=/Share2/home/zhangqf7/yuhan/6_Docking/20220225_datasets_PDB/$1
DRUG=/Share2/home/zhangqf7/yuhan/6_Docking/20220225_datasets_PDB/Drug_pdb

file=$2
filename=$(echo $file | cut -d . -f1)
final_data=${OUT}/$3


cd ${OUT}
mkdir ${final_data}
chimera --nogui --script "${TOOLS}/dockprep.py ${DATA}/${file} ${OUT}"
cd ${OUT}/$filename
sh ${TOOLS}/box.sh
c_x=$(sed -n '2p' 2_Rec_box.pdb | awk '{print $6}')
c_y=$(sed -n '2p' 2_Rec_box.pdb | awk '{print $7}')
c_z=$(sed -n '2p' 2_Rec_box.pdb | awk '{print $8}')
s_x=$(sed -n '3p' 2_Rec_box.pdb | awk '{print $6}')
s_y=$(sed -n '3p' 2_Rec_box.pdb | awk '{print $7}')
s_z=$(sed -n '3p' 2_Rec_box.pdb | awk '{print $8}')
/Share2/home/zhangqf7/yuhan/6_Docking/autodock/results_pdb/programs/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ${Utilities}/prepare_receptor4.py -r 1_Rec_noH.pdb -o 1_Rec_noH.pdbqt -U waters


for ligand in `ls ${DRUG} --color=never | sort -V`
do
ligname=$(echo $ligand | cut -d . -f1)

mkdir ${OUT}/$filename/Ligands/${ligand}
cd ${OUT}/$filename/Ligands/${ligand}

chimera --nogui --script "${TOOLS}/ligprep.py ${DRUG}/${ligand} ${OUT}"

timeout 60m /Share2/home/zhangqf7/yuhan/6_Docking/autodock/results/programs/mgltools_x86_64Linux2_1.5.6/bin/pythonsh ${TOOLS}/prepare_ligand4.py -l 1_Lig.pdb -o 1_Lig.pdbqt

timeout 180m bash ${TOOLS}/docking.sh ../../1_Rec_noH.pdbqt 1_Lig.pdbqt $c_x $c_y $c_z $s_x $s_y $s_z ${TOOLS}/rmH.py > 3_predictions_out.txt

cp 3_predictions_scored_1.pdb 3_predictions_scored_best.pdb

weight=$(obprop 3_predictions_scored_best.pdb | grep ^mol_weight | awk '{print $2}')
#encode=$(obprop 3_predictions_scored_best.pdb | grep ^canonical_SMILES | awk '{print $2}')

#binding affinity
grep + -A 50 3_predictions_out.txt | awk '$1==1{print $0}' | awk '{print $2}' | awk -v rec_name="${filename}" -v lig_name="${ligand}" -v lig_weight="${weight}"  '{print rec_name "\t" lig_name "\t" lig_weight "\t" $1 }' > 3_predictions_score_best.txt

chimera --nogui --nostatus --script "${TOOLS}/combine.py ${OUT}/$filename/1_Rec_noH.pdb 3_predictions_scored_best.pdb ${TOOLS}" #4_results_final.pdb
chimera --nogui --nostatus --script "${TOOLS}/rmH.py 4_results_final.pdb 4_results_final_noH.pdb"
chimera --nogui --nostatus --script "${TOOLS}/combine2.py ${OUT}/$filename/1_Rec_noH.pdb 3_predictions_scored_best.pdb" #4_results_tmp.pdb

python ${TOOLS}/merge.py -i 4_results_final_noH.pdb -r 4_results_tmp.pdb -o 5_results.pdb

score=$(grep + -A 50 3_predictions_out.txt | awk '$1==1{print $0}' | awk '{print $2}' | bc -l | xargs printf "%.1f")

cp 5_results.pdb ${final_data}/${filename}_${ligname}_${score}.pdb

cat ${OUT}/$filename/Ligands/${ligand}/3_predictions_score_best.txt >> ${OUT}/$filename/3_predictions_score_best_all.txt

done


