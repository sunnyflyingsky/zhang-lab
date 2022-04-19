tranFa=../ref/hg38_transcriptome.fa

perl ../ref/seqMerge.pl 2 ../Step1_dynamic/window-anno.bed ${tranFa} > fg.fa

nohup fimo --oc fimo_out --thresh 0.001 ../ref/Homo_sapiens.meme fg.fa &

for i in `seq 1 100`
do
	shuffleBed -i ../Step1_dynamic/window-anno.bed -chrom -g ../ref/transcript-length.txt -seed ${i} > window-anno-shuffle.bed
	rm ../ref/hg38_transcriptome.fa.fai
	bedtools getfasta -fi ${tranFa} -bed window-anno-shuffle.bed -fo fg_shuffle.fa
	fimo --norc --oc res/fimo_out_shuffle_${i} --thresh 0.001 ../ref/Homo_sapiens.meme fg_shuffle.fa
	wait;
	cp res/fimo_out_shuffle_${i}/fimo.tsv res/fimo_${i}.tsv
	rm -rf res/fimo_out_shuffle_${i}
	wait;
done

cat res/*.tsv > fimo2.txt
cp fimo_out/fimo.tsv fimo1.txt

#perl ../Step2_generate/seqMerge.pl 2 window-anno-shuffle.bed ${tranFa} > fg_shuffle.fa
#nohup fimo --oc fimo_out_shuffle --thresh 0.001 Collapsed.used.meme fg_shuffle.fa &

