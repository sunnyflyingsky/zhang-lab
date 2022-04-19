cd /data/yuhan/dssr_details
for file in `ls /data/yuhan/dssr_details/PDB_RNA_all | sort -n`
do
filename=$(echo $file | cut -d . -f1)
echo ${filename} ;
/home/yuhan/dssr/x3dna-dssr --input=/data/yuhan/dssr_details/PDB_RNA_all/${filename}.pdb --output=/data/yuhan/dssr_details/results/${filename}.out
num=$(grep structural /data/yuhan/dssr_details/results/${filename}.out | awk '{print $6}')
grep -A ${num} linkage /data/yuhan/dssr_details/results/${filename}.out | sed -e '1d' > /data/yuhan/dssr_details/${filename}.txt
done


#bulge
#hairpin-loop
#internal-loop
#junction-loop
#ss-non-loop
