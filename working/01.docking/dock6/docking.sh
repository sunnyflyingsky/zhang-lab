infile="/Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA"
pythonScriptPath="/Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/scripts"

find ${infile}/data_set -name *.pdb > temp.txt
while read line
do
    pdbFile=$(echo $line)
	name=${pdbFile:0-8:4}
	chimera --nogui --script "ChimeraPrep.py -i "${pdbFile}" -o "${infile}"/temp_set"
	wait;
	eval $(awk '
	{if($1=="x_lenghth:"){l = $2} 
	else if($1=="y_width:"){w = $2} 
	else if($1=="z_height:"){h = $2} 
	else if($1=="center:"){cx = $2;cy = $3;cz = $4}}
	END{printf("l=%s;w=%s;h=%s;cx=%s;cy=%s;cz=%s",l,w,h,cx,cy,cz)}' ${infile}/temp_set/windows.txt)
	wait;
	cd temp_set/site_file
	sphgen INSPH
	showsphere < sphgen_cluster.in
	sphere_selector receptor.sph ${infile}/temp_set/struct_file/ligand.mol2 10.0
	showsphere < selected_spheres.in
	rm OUTSPH
	rm receptor.sph
	wait;
	cd ../generation_grid
	showbox < box.in
	grid -i grid.in -o grid.out
	wait;
	cd ../lig_sampling
	dock6 -i rigid.in -o rigid.out
	wait;
	dock6 -i anchor_and_grow.in -o anchor_and_grow.out
	wait;
	cd ../..
	cp temp_set/lig_sampling/anchor_and_grow_ranked.mol2 res/${name}/${name}_anchor_and_grow_ranked.mol2
	cp temp_set/lig_sampling/rigid_ranked.mol2 res/${name}/${name}_rigid_ranked.mol2
	cp temp_set/lig_sampling/anchor_and_grow.out res/${name}/${name}_anchor_and_grow.out
	cp temp_set/lig_sampling/rigid.out res/${name}/${name}_rigid.out
	rm temp_set/lig_sampling/rigid.out
	rm temp_set/lig_sampling/anchor_and_grow.out

done < ${infile}/temp.txt
	