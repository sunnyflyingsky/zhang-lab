infile="/Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA"
pythonScriptPath="/Share2/home/zhangqf7/yuhan/rotation_student/zhangjiasheng/Docking/RNA/scripts"
pythonshPath="/Share2/home/zhangqf7/yuhan/6_Docking/autodock/results/programs/mgltools_x86_64Linux2_1.5.6/bin"

find ${infile}/data_set -name *.pdb > temp.txt
while read line
do
    pdbFile=$(echo $line)
	name=${pdbFile:0-8:4}
	chimera --nogui --script "ChimeraPrep.py -i "${pdbFile}" -o "${infile}"/temp_set"
	wait;
	cd temp_set
	sphgen INSPH
	sphere_selector receptor.sph ${infile}/temp_set/struct_file/ligand.mol2 20.0
	showsphere < selected_spheres.in
	rm OUTSPH
	rm receptor.sph
	showbox < box.in
	cd ..
	wait;

	eval $(awk '
	{if($1=="x_lenghth:"){l = $2} 
	else if($1=="y_width:"){w = $2} 
	else if($1=="z_height:"){h = $2} 
	else if($1=="center:"){cx = $2;cy = $3;cz = $4}}
	END{printf("l=%s;w=%s;h=%s;cx=%s;cy=%s;cz=%s",l,w,h,cx,cy,cz)}' ${infile}/temp_set/windows.txt)
	${pythonshPath}/pythonsh ${pythonScriptPath}/prepare_receptor4.py -r ${infile}/temp_set/${name}_receptor.pdb -o ${infile}/temp_set/${name}_receptor.pdbqt -A checkhydrogens -U waters
	${pythonshPath}/pythonsh ${pythonScriptPath}/prepare_ligand4.py -l ${infile}/temp_set/${name}_ligand.pdb -o ${infile}/temp_set/${name}_ligand.pdbqt
	wait;
	vina --receptor ${infile}/temp_set/${name}_receptor.pdbqt --ligand ${infile}/temp_set/${name}_ligand.pdbqt \
		--center_x ${cx} --center_y ${cy} --center_z ${cz} --size_x ${l} --size_y ${w} --size_z ${h} \
		--out ${infile}/res/${name}.pdbqt \
		--exhaustiveness 10 --num_modes 10 --energy_range 1.36
	wait;
	mkdir ${infile}/res/${name}
	babel -i pdbqt ${infile}/res/${name}.pdbqt -o pdb ${infile}/res/${name}/${name}.pdb
	babel -i pdbqt ${infile}/temp_set/${name}_ligand.pdbqt -o pdb ${infile}/res/${name}/${name}_ligand.pdb
	chimera --nogui --script "cal_rmsd.py -i "${infile}"/temp_set/"${name}".pdb -o "${infile}"/res/"${name}
	wait;
done < ${infile}/temp.txt
	