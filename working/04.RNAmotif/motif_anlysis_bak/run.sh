## run this file 

cd Step2_generate
bash run.sh
wait;
cd ../Step3_fimo_shuffle
bash shuffle.sh
cp fimo1.txt ../Step4_enrich/fimo1.txt
cp fimo2.txt ../Step4_enrich/fimo2.txt
wait;
cd ../Step4_enrich
python enrich.py
wait;
cd ../Step5_plot
nohup python trans_meme.py &
nohup python F2.motif_scatter_heatmap.py &
wait;
