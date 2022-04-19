运行脚本文件，dynamicRegion.py,
/data/yuhan/zhangjiasheng/motifEnrichment/dynamic_region
运行路径下存在文件目录下所有脚本文件

    <Requires>
    -i parameter.in

    you should input the input information in parameter.in, \n
    including: 
    inputFolder, all of your input file should be in this folder path, including eg1,eg2,egF,cg1,cg2,cgF.
    refFolder, all of your annotation file should be in this folder path, including genomeCoor and transcriptome file.
    ref_genomeCoor_file, genomeCoor file.
    ref_transcriptome_file, transcriptome file.
    outputFolder, we will output the results to this path. and if you start from step>1(if you first use, please set the start_step to 1), please put the file in this folder, including dynamic_region_demo, filtered, merged, anno, etc.  
    inputEgFile1, experimental icshape score, group1, for calculate difference.
    inputEgFile2, experimental icshape score, group2, for calculate difference.
    inputCgFile1, control icshape score, group1, for calculate difference.
    inputCgFile2, control icshape score, group2, for calculate difference.
    inputEgFileFinal, experimental icshape score for calculate dynamic region
    inputCgFileFinal, control icshape score for calculate dynamic region

    [Options]
    dynamicWindows, window size for dynamic search
    dynamicSteps, window step for dynamic search
    slidingWindows, slide window size for cutting the merged region to get same size regions 
    slidingSteps, slide window step for cutting the merged region to get same size regions
    filter_cutoff, filter cutoff for P-value
    None_cutoff, cutoff for label the region with much None score siganl. the larger it is, the more region will be labeled unknown.  
    struct_cutoff, cutoff for label the region to more single-structure or less single-structure. the larger it is, the more region will be labeled unknown. 
    start_tep, which step you want to start. if you first use, please set it to 1. it should be a int number from 1 to 6.  
                1 means from shape file to get all the results file
                2 means from dynamic_region_demo file to get onther results file
                3 means from dynamic_region_filtered file to get following results file
                4 means from dynamic_region_filtered_merged file to get following results file
                5 means from dynamic_region_filtered_merged_slided file to get following results file
                6 means from dynamic_region_filtered_merged_slided_anno file to get final window-anno.bed file
    dynamic_positive_num, cutoff for select dynamic windows, which have more sites(score>indiff) than this num. 
    dynamic_valid_num, cutoff for select dynamic windows, which have more valid score than this num.
    UTR3_or_ALLREGION, if you want to get UTR3 region. True for UTR3, False for all region.

	ex: python dynamicRegion.py -i parameter.in
	
	[default]
	inputFolder             /data/yuhan/zhangjiasheng/motifEnrichment/dynamic_region_analysis_NES_ERM
	refFolder               /data/yuhan/zhangjiasheng/motifEnrichment/dynamic_region/ref
	outputFolder            /data/yuhan/zhangjiasheng/motifEnrichment/dynamic_region/res
	inputEgFile1            NES_R1_FPKM.shape
	inputEgFile2            NES_R2_FPKM.shape
	inputCgFile1            ERM_R1_FPKM.shape
	inputCgFile2            ERM_R2_FPKM.shape
	inputEgFileFinal        NES_FPKM.shape
	inputCgFileFinal        ERM_FPKM.shape
	ref_genomeCoor_file     hg38.genomeCoor.bed
	ref_transcriptome_file  hg38_transcriptome.fa
	dynamicWindows          101
	dynamicSteps            1
	slidingWindows          30
	slidingSteps            10
	filter_cutoff           0.0001
	None_cutoff             0.6
	struct_cutoff           0
	start_step              1
	dynamic_positive_num    5
	dynamic_valid_num       60
	UTR3_or_ALLREGION       True

输入文件 parameter.in
运行环境 conda activate icshape

parameter文件的参数解释
	inputFolder, input的数据的路径，一般包括六个文件，两个实验组的shape，两个对照组的shape，以及合并后的实验组和对照组shape文件。
    refFolder, 参考基因组的数据路径，路径下包括两个文件genomeCoor以及transcriptome的文件
    ref_genomeCoor_file, genomeCoor文件的名字
    ref_transcriptome_file, transcriptome文件的名字
    outputFolder, 输出文件路径，我们会将所有逐步输出的文件输出在此路径下，且输出文件的name保持一致。如果你在start_step中选取了从后面的步骤起始，请保证每一个step number对应的起始文件在输出路径下，对应文件我们将在start_step中详细介绍
			输出文件包括：
			results.txt
			results_filter.txt
			results_merge_filter.txt
			results_merge_filter_window.txt
			results_merge_filter_window_anno.txt
			results_merge_filter_window_anno_utr3.txt
			window-anno.bed
	inputEgFile1, 实验组shape数据
    inputEgFile2, 实验组shape数据
    inputCgFile1, 对照组shape数据
    inputCgFile2, 对照组shape数据
    inputEgFileFinal, 合并后的实验组shape数据
    inputCgFileFinal, 合并后的对照组shape数据

    [Options]
    dynamicWindows, 进行dynamic region搜索时的滑动窗口大小(bp)
    dynamicSteps, 进行dynamic region搜索时的滑动窗口滑动距离(bp)
    slidingWindows, 切割窗口的大小，切割merge后的不等长的dynamic region到相同大小
    slidingSteps, 切割窗口的滑动距离，每次切割到下一次切割起始位置的距离
    filter_cutoff, 过滤dynamic region的cutoff值，对P-value进行筛选
    None_cutoff, 用于将一段含有过多None值的序列标记为Unknown的阈值，当None值多于这个阈值时，我们会将这段序列的structure倾向性标记为unknown 
    struct_cutoff, 用于评判一段序列的单链倾向性和双链倾向性。
				差异性shape的差值正值数量-负值数量多于这个阈值时，我们认为这段序列结构倾向性为more single structure，标记为“more” 
				差异性shape的差值负值数量-正值数量多于这个阈值时，我们认为这段序列结构倾向性为less single structure，标记为“less”
				当处于中间时（即不满足上述两种条件时），我们认为这段序列结构不具备倾向性，标记为“unknown”
    start_tep, 起始步骤，表示您希望从哪一步起始。我们用输出文件来分隔程序，一种有六种起始选择，1,2,3,4,5,6. 小于1和大于6的值不被允许，我们会默认重置为1。
				由于各个输出文件结构的问题，请务必保证文件结构正确，如果是第一次运行，请设置为1，如果你清楚这些文件，请按照自己的需求设置。
                1 代表从最开始运行，包含从dynamic region的搜索以及往后的步骤。输入文件为在inputFolder路径下的六个shape文件。输出文件到outputFolder路径下，一共有7个。
                2 代表从筛选region以及往后的步骤。输入文件为在outputFolder路径下的results.txt文件。输出文件到outputFolder路径下，一共有6个。
                3 代表从合并region以及往后的步骤。输入文件为在outputFolder路径下的results_filter.txt文件。输出文件到outputFolder路径下，一共有5个。
                4 代表从切割region以及往后的步骤。输入文件为在outputFolder路径下的results_merge_filter.txt文件。输出文件到outputFolder路径下，一共有4个。
                5 代表从注释region以及往后的步骤。输入文件为在outputFolder路径下的results_merge_filter_window.txt文件。输出文件到outputFolder路径下，一共有3个。
                6 代表从截取UTR3以及往后的步骤。输入文件为在outputFolder路径下的results_merge_filter_window_anno.txt文件。输出文件到outputFolder路径下，一共有2个。
    dynamic_positive_num, 用于进行dynamic region检索时进行筛选，仅当差异分数的绝对值大于组内差异的位点会被计数，计数超过这个阈值的region会被输出保留
    dynamic_valid_num, 用于进行dynamic region检索时进行筛选，仅当非None的shape值计数超过这个阈值时，这个region才会被认为有效。
    UTR3_or_ALLREGION, 是否输出UTR3的dynamic region文件。
				True，此时输出七个文件，包含results_merge_filter_window_anno_utr3.txt，然后window-anno.bed由UTR3的文件生成
				False，此时输出六个文件，不包含results_merge_filter_window_anno_utr3.txt。window-anno.bed由results_merge_filter_window_anno.txt生成，包括所有区域的dynamic region。



