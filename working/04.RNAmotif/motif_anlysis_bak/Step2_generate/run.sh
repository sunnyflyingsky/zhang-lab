
	script=../ref/seqMerge.pl 

	#########区分utr5 cds utr3
	grep "utr5" ../Step1_dynamic/window-anno.bed> window-anno_utr5.bed
	grep "cds" ../Step1_dynamic/window-anno.bed> window-anno_cds.bed
	grep "utr3" ../Step1_dynamic/window-anno.bed> window-anno_utr3.bed
	
	#######取序列做motif分析
	script=../ref/seqMerge.pl 
	tranFa=../ref/hg38_transcriptome.fa
	
	perl $script 2 ../Step1_dynamic/window-anno.bed $tranFa  > fg.fa
	perl $script 2 window-anno_utr5.bed $tranFa  > fg_utr5.fa
	perl $script 2 window-anno_cds.bed $tranFa  > fg_cds.fa
	perl $script 2 window-anno_utr3.bed $tranFa  > fg_utr3.fa

