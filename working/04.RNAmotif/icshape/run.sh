refPath=
inPath=

##################################################
MOUSE_ABUNDANT_RNA=/Share2/home/zhangqf7/yuhan/0_genome/bowtie2/mouse_rRNA_tRNA_mtRNA/mouse_rRNA_tRNA_mtRNA  ## Abundant RNAs
GENOME=/Share2/home/zhangqf7/yuhan/0_genome/STAR/mm10_Gencode
ANNOTATION=/Share2/home/zhangqf7/yuhan/0_annotation/Gencode/mm10.genomeCoor.bed

HUMAN_ABUNDANT_RNA=/Share2/home/zhangqf7/yuhan/0_genome/bowtie2/human_rRNA_tRNA_mtRNA/human_rRNA_tRNA_mtRNA
GENOME=/Share2/home/zhangqf7/yuhan/0_genome/STAR/hg38_Gencode
ANNOTATION=/Share2/home/zhangqf7/yuhan/0_annotation/Gencode/hg38.genomeCoor.bed

## bsub -o output.txt -e err.log -q Z-ZQF -m node531 -n 1 "xxx"

###################################################

## 创建软连接和硬链接
ln -s /Share2/home/zhangqf7/yuhan/1_xiwen/20211213_TAM/0.rawdata/cKO_R2.fq 01.rawdata/cKO_R2.fq
ln /Share2/home/zhangqf7/yuhan/1_xiwen/20211213_TAM/0.rawdata/cKO_R2.fq 01.rawdata/cKO_R2.fq

##构建基因组注释和索引
bowtie2-build --quiet /Share2/home/zhangqf7/yuhan/0_genome/bowtie2/human_rRNA_tRNA_mtRNA/human_rRNA_tRNA_mtRNA.fa ref/rRNA_index/rRNA
falen /Share2/home/zhangqf7/yuhan/0_genome/bowtie2/human_rRNA_tRNA_mtRNA/human_rRNA_tRNA_mtRNA.fa > ref/rRNA_index/rRNA.len

icSHAPE-pipe starbuild -i ref/hg38.fa -o ref/index/ --gtf ref/hg38.gtf -p 20
icSHAPE-pipe parseGTF -g ref/hg38.gtf -o ref/GTF/Anno -s ensembl --genome ref/hg38.fa
fetchSmallRNA.py ref/GTF/Anno_transcriptome.fa ref/smallRNA/smallRNA.fa 200
bowtie2-build --quiet ref/smallRNA/smallRNA.fa ref/smallRNA/smallRNA
falen ref/smallRNA/smallRNA.fa > ref/smallRNA/smallRNA.len

##去重相同序列的reads，避免PCR duplicates
似乎在这一步丢失了barcode信息，这个信息是否重要
icSHAPE-pipe readcollapse -U 01.rawdata/cKO_R2.fq -o 02.readCollapse/cKO_R2.fq --simplify

##Trim 3'adaptor与5' random barcord
不清楚icSHAPE与smartSHAPE是否存在区别
fastqc结果表明，trim后数据质量提升，但是依然存在一些小瑕疵，目前没有找到原因(其中icSHAPE数据表现较好，但是smartSHAPE数据出现duplicate异常)
为什么文献中提及的是13nt。师兄说的是10nt
我所进行trim的adaptor去除了poly尾端，最终结果会比师兄的略少(03.trim/.fq)
icSHAPE-pipe trim -i 02.readCollapse/cKO_R2.fq -o 03.trim/cKO_R2.fq -l 10 -a ref/adaptor.fa -p 20 -m 25

##去除map到rRNA和smallRNA上的reads
icSHAPE-pipe cleanFq -i 03.trim/G3BP1_NBP_R1.fq -o 04.rem_rRNA/G3BP1_NBP_R1.fq -x ref/rRNA_index/rRNA -p 20 --mode Local --sam 04.rem_rRNA/G3BP1_NBP_R1.sam
icSHAPE-pipe cleanFq -i 04.rem_rRNA/G3BP1_NBP_R1.fq -o 04.rem_smallRNA/G3BP1_NBP_R1.fq -x ref/smallRNA/smallRNA -p 20 --mode Local --sam 04.rem_samllRNA/G3BP1_NBP_R1.sam

##mapping
icSHAPE-pipe mapGenome -i 04.rem_smallRNA/G3BP1_DBP_R1.fq -o 05.mapGenome/G3BP1_DBP_R1 -x ref/index -p 20 --noMut5


##计算FPKM
为什么只需要对背景计算
DBP
icSHAPE-pipe calcFPKM -i 05.mapGenome/G3BP1_DBP_R2.sorted.bam -o 06.calcFPKM/G3BP1_DBP_R2 -G ref/hg38.gtf -p 20

##计算score
icSHAPE-pipe sam2tab -in 05.mapGenome/G3BP1_NBP_R1.sorted.bam -out 07.sam2tab/G3BP1_NBP_R1.tab
icSHAPE-pipe sam2tab -in 04.rem_rRNA/G3BP1_NBP_R1.sam -out 07.sam2tab/rRNA_G3BP1_NBP_R1.tab
icSHAPE-pipe sam2tab -in 04.rem_smallRNA/G3BP1_NBP_R1.sam -out 07.sam2tab/smallRNA_G3BP1_NBP_R1.tab


##icSHAPE score计算
为什么是ATCG，而不是AUCG
N1,N2有次序的吗
small_RNA没有结果，为什么
icSHAPE-pipe calcSHAPE \
	-D 07.sam2tab/G3BP1_DBP_R1.tab,07.sam2tab/G3BP1_DBP_R2.tab \
	-N 07.sam2tab/G3BP1_NBP_R1.tab,07.sam2tab/G3BP1_NBP_R2.tab \
	-size ref/index/chrNameLength.txt \
	-ijf ref/index/sjdbList.fromGTF.out.tab \
	-genome ref/hg38.fa \
	-bases A,T,C,G \
	-out 08.calcGenomeSHAPE/shape.gTab

icSHAPE-pipe calcSHAPE \
	-D 07.sam2tab/rRNA_G3BP1_DBP_R1.tab,07.sam2tab/rRNA_G3BP1_DBP_R2.tab \
	-N 07.sam2tab/rRNA_G3BP1_NBP_R1.tab,07.sam2tab/rRNA_G3BP1_NBP_R2.tab \
	-size ref/rRNA_index/rRNA.len \
	-genome ref/human_rRNA_tRNA_mtRNA.fa \
	-bases A,T,C,G \
	-out 08.calcGenomeSHAPE/rRNA_shape.gTab \
	-non-sliding

icSHAPE-pipe calcSHAPE \
	-D 07.sam2tab/smallRNA_G3BP1_DBP_R1.tab,07.sam2tab/smallRNA_G3BP1_DBP_R2.tab \
	-N 07.sam2tab/smallRNA_G3BP1_NBP_R1.tab,07.sam2tab/smallRNA_G3BP1_NBP_R2.tab \
	-size ref/smallRNA/smallRNA.len \
	-genome ref/human_rRNA_tRNA_mtRNA.fa \
	-bases A,T,C,G \
	-out 08.calcGenomeSHAPE/smallRNA_shape.gTab \
	-non-sliding

## 基因组icSHAPE score转成转录本icSHAPE score
-r参数的意义
最后的到的final.shape结果很少，为什么，是注释文件有误还是中间某一步出错了
icSHAPE-pipe genSHAPEToTransSHAPE \
	-i 08.calcGenomeSHAPE/shape.gTab \
	-o 09.transSHAPE/final.shape \
	-g ref/GTF/Anno.genomeCoor.bed \
	-p 20 \
	-c 100 \
	-T 0.1 \
	-m 0.8 \ 
	-r 06.calcFPKM/G3BP1_DBP_R1/isoforms.fpkm_tracking,06.calcFPKM/G3BP1_DBP_R2/isoforms.fpkm_tracking

icSHAPE-pipe genSHAPEToTransSHAPE \
	-i 08.calcGenomeSHAPE/rRNA_shape.gTab \
	-o 09.transSHAPE/rRNA_final.shape \
	-s ref/rRNA_index/rRNA.len \
	-p 1 \
	--app

## 可视化
icSHAPE-pipe genSHAPEToBedGraph -i 08.calcGenomeSHAPE/shape.gTab -o 10.bedGraph/ -c 200


#####################################################################################################

icSHAPE-pipe calcSHAPE \
	-D ${OUT}/6.sam2tab/G3BP1_ADBH_R1.fq.sorted.bam.tab \
	-N ${OUT}/6.sam2tab/G3BP1_ANBH_R1.fq.sorted.bam.tab  \
	-size ${GENOME}/chrNameLength.txt \
	-ijf ${GENOME}/sjdbList.fromGTF.out.tab \
	-genome ${GENOME}/hg38.fa -bases A,T,C,G \
	-omc 100 \
	-out ${OUT}/7.calcGenomeSHAPE/shape_A_R1.gTab

icSHAPE-pipe genSHAPEToTransSHAPE \
	-i ${OUT}/7.calcGenomeSHAPE/shape_A_R1.gTab \
	-o ${OUT}/8.transSHAPE/final_A_R1.shape \
	-g ${ANNOTATION} \
	-p 5 \
	-c 100 \
	-T 0.01 \
	-m 0.8
##不采用-r参数去除异构型，得到shape_more文件










