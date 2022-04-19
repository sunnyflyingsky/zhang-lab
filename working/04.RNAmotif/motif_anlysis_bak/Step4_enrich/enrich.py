import subprocess
import random
from nested_dict import nested_dict
import pandas as pd
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests

def read_fimo_txt(txt=None):
	fimo_dict = nested_dict(1, list)
	n = 0
	with open(txt, 'r') as TXT:
		for line in TXT:
			line = line.strip()
			if line.startswith('#') or not line:
				continue
			n += 1
			arr = line.split('\t')
			fimo_dict[arr[0]].append(arr[2])
			#fimo_dict[arr[0]].append(arr[2])  ## 2021/12/12 zjs True solving count problems
	for i,j in fimo_dict.items(): # one RBP may bind multiple site of a region
		#fimo_dict[i] = list(j)  ##2021/12/09 zjs count problems
		fimo_dict[i] = list(set(j))  
	print("read: %s"%(txt))
	print("line: %s"%(n))
	print("RBPs: %s"%(len(fimo_dict)))
	print({i:len(j) for i,j in fimo_dict.items()})
	print()
	return fimo_dict.to_dict()

def read_RBP_name():
	txt = 'RNAmotif/motif_anlysis/Step4_enrich/RBP_human.txt'
	df = pd.read_csv(txt, header=None, sep='\t')
	return {i:j for i,j in zip(df[1], df[2])}

def random_select():
	pass

def fimo_txt_compare_enrich(fimo_txt1=None, fimo_txt2=None, savefn=None):
	#if fimo_txt1 is None:
	#	fimo_txt1 = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/new_mergepeaks_d10/separate/way2345.sorted.merge.annoused.utr3/fimo.tsv'
	#if fimo_txt2 is None:
	#	fimo_txt2 = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/new_mergepeaks_d10/separate/shuffle/way2345.sorted.merge.annoused.utr3.s12345/fimo.tsv'
	#if savefn is None:
	#	savefn = fimo_txt2.replace('.tsv', '.enrich.txt')

	# bed1 = fimo_txt1.replace('/fimo.tsv', '.bed')
	# bed2 = fimo_txt2.replace('/fimo.tsv', '.bed')

	bed1 = fimo_txt1
	bed2 = fimo_txt2

	df_bed1 = pd.read_csv(bed1, header=None, sep='\t')
	df_bed2 = pd.read_csv(bed2, header=None, sep='\t')
	fimo_dict1 = read_fimo_txt(fimo_txt1)
	fimo_dict2 = read_fimo_txt(fimo_txt2)

	RBP_common = set(fimo_dict1.keys()) & set(fimo_dict2.keys())
	print(RBP_common)
	RBP_dict = nested_dict()

	for RBP in RBP_common:
		RBP_dict[RBP]['RBP'] = RBP
		fimo1_motif =  len(fimo_dict1[RBP])
		RBP_dict[RBP]['fimo1_motif'] = fimo1_motif
		fimo1_not_motif = df_bed1.shape[0] - len(fimo_dict1[RBP])
		RBP_dict[RBP]['fimo1_not_motif'] = fimo1_not_motif
		fimo2_motif = len(fimo_dict2[RBP])
		RBP_dict[RBP]['fimo2_motif'] = fimo2_motif
		fimo2_not_motif = df_bed2.shape[0] - len(fimo_dict2[RBP])
		RBP_dict[RBP]['fimo2_not_motif'] = fimo2_not_motif
		RBP_dict[RBP]['odd1'] = fimo1_motif / float(fimo1_not_motif)
		RBP_dict[RBP]['odd2'] = fimo2_motif / float(fimo2_not_motif)
		#print(RBP, [[fimo1_motif, fimo1_not_motif], [fimo2_motif, fimo2_not_motif]])
		oddsratio, pvalue = stats.fisher_exact([[fimo1_motif, fimo1_not_motif], [fimo2_motif, fimo2_not_motif]])
		#print(pvalue)
		RBP_dict[RBP]['odd'] = oddsratio
		RBP_dict[RBP]['pvalue'] = pvalue

	print(RBP_dict)
	RBP_df = pd.DataFrame.from_dict(RBP_dict, orient='index')
	print(RBP_df.head())
	_, pvalue_adj, _, _ = multipletests(list(RBP_df['pvalue']), method='bonferroni')
	#print(pvalue_adj)
	#print(RBP_df.shape, len(pvalue_adj))
	RBP_df['pvalue_adj'] = pvalue_adj
	print(RBP_df.head())

	cols_order = ['RBP', 'RBP_name_simple','RBP_name_complex', 'fimo1_motif', 'fimo1_not_motif', 'fimo2_motif', 'fimo2_not_motif', 'odd1', 'odd2', 'odd', 'pvalue', 'pvalue_adj']
	RBP_name_dict = read_RBP_name()
	#RBP_df['RBP_name'] = [RBP_name_dict[i] if RBP_name_dict.has_key(i) else i for i in RBP_df['RBP']]
	#RBP_df_new = RBP_df.copy()
	#RBP_df['RBP_name'] = [RBP_name_dict[i] if i in RBP_name_dict else i for i in RBP_df['RBP']]
	RBP_df['RBP_name_complex'] = [RBP_name_dict[i] if i in RBP_name_dict else i for i in RBP_df['RBP']]
	RBP_df['RBP_name_simple'] = [RBP_name_dict[i].split('_')[0] if i in RBP_name_dict else i for i in RBP_df['RBP']]
	RBP_df[cols_order].to_csv(savefn, header=True, sep='\t', index=False)
	#RBP_df_new[cols_order].to_csv(savefn.replace('/fimo.enrich.txt','/fimo_new.enrich.txt'), header=True, sep='\t', index=False)


#read_path = 'RNAmotif/motif_anlysis/Step4_enrich/'
read_path = '../Step4_enrich/'
fimo_txt1 = read_path+'fimo1.txt'
fimo_txt2 = read_path+'fimo2.txt'
savefn = read_path+'fimo.enrich.txt'
#fimo_txt1 = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/new_mergepeaks_d10/separate/way2345.sorted.merge.annoused.utr3/fimo.tsv'
#fimo_txt2 = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_merge_region/005_005_new/abs/new_mergepeaks_d10/separate/shuffle/way2345.sorted.merge.annoused.utr3.s12345/fimo.tsv'
fimo_txt_compare_enrich(fimo_txt1, fimo_txt2, savefn)
