#rna2meme 
function FIMO(){
	cd $1
	fimo=~/bin/fimo
	fg=fg.fa
	fg_utr5=fg_utr5.fa
	fg_cds=fg_cds.fa
	fg_utr3=fg_utr3.fa

	motif=/Share/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Working/04.RNAmotif/Step3_fimo_shuffle/Homo_sapiens.meme
	
	nohup $fimo --oc fimo_out-all  --thresh 0.001 $motif  $fg &
	nohup $fimo --oc fimo_out-utr5  --thresh 0.001 $motif  $fg_utr5 &
	nohup $fimo --oc fimo_out-cds  --thresh 0.001 $motif  $fg_cds &
	nohup $fimo --oc fimo_out-utr3  --thresh 0.001 $motif  $fg_utr3 &

}

for j in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_sphere" "sphere_shield"
do
		dir=$1/up/$j
		echo $dir
		FIMO  $dir &
done


for j in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_sphere" "sphere_shield"
do
		dir=$1/down/$j
		echo $dir
		FIMO  $dir &
done


for j in "egg_cell1"  "cell1_cell4" "cell4_cell64" "cell64_sphere" "sphere_shield"
do
		dir=$1/abs/$j
		echo $dir
		FIMO  $dir &
done




	




	
	
	
	
	
	
