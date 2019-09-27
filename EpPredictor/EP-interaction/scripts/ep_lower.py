import sys
sys.path.append('../../')


import numpy as np
from Bio import SeqIO
import pandas as pd
import argparse


from core import utils

use_cols = ['chrom-Enh','chromStart','chromEnd','TSS','label','id']


def _seqmat2vec(seq_mat):
	"""
	@seq_mat: [[seq],[seq],[seq]], numpy array
	"""

	_, cols = seq_mat.shape
	init_mat = np.zeros((4,cols))

	for col in range(cols):
		col_2_list = seq_mat[:,col].tolist()
		count_A = col_2_list.count('A')
		count_T = col_2_list.count('T')
		count_C = col_2_list.count('C')
		count_G = col_2_list.count('G')
		tmp = [count_A,count_T,count_C,count_G]
		init_mat[:,col] = tmp
	init_mat = init_mat.reshape(1,-1)

	return init_mat[0]

def _sequence_region2mat(sequence_string, feature_cols):
	total_mat = []
	init_lement_idx = 0
	while(init_lement_idx<=len(sequence_string)):
		tmp = list(sequence_string[init_lement_idx:init_lement_idx+feature_cols])
		total_mat.append(tmp)
		init_lement_idx+=feature_cols
	total_mat = pd.DataFrame(total_mat)
	total_mat = total_mat.fillna('F')
	return total_mat.values

def generate_feature(filepath,feature_cols,required_info_df,chr_name):

	feature_df = []

	record_seq = None
	for record in SeqIO.parse(filepath, "fasta"):
		record_seq = record.seq

	for each_std_key in required_info_df.values:
		info = []
		
		enahncer_start = each_std_key[1]
		enhancer_end = each_std_key[2]
		tss = each_std_key[3]
		label = each_std_key[4]
		id = each_std_key[-1]

		sequence_enhancer = str(record_seq[enahncer_start:enhancer_end]).upper()
		if(enhancer_end>tss):
			enhancer_end,tss = tss, enhancer_end
		sequence_window = str(record_seq[enhancer_end:tss]).upper()
		sequence_patch = str(record_seq[enahncer_start:tss]).upper()
		sequence_promoter2k = str(record_seq[tss-2000:tss+2000]).upper()

		info = [chr_name, enahncer_start, enhancer_end,tss, label,id]

		sequence_window = _sequence_region2mat(sequence_window, feature_cols)
		sequence_enhancer = _sequence_region2mat(sequence_enhancer, feature_cols)
		sequence_patch = _sequence_region2mat(sequence_patch, feature_cols)
		sequence_promoter2k = _sequence_region2mat(sequence_promoter2k, feature_cols)
		

		sequence_window = _seqmat2vec(sequence_window)
		sequence_enhancer = _seqmat2vec(sequence_enhancer)
		sequence_patch = _seqmat2vec(sequence_patch)
		sequence_promoter2k = _seqmat2vec(sequence_promoter2k)
		info_all = np.r_[info,sequence_enhancer, sequence_window,sequence_patch,sequence_promoter2k]
		feature_df.append(info_all)

	feature_df = pd.DataFrame(feature_df)

	return feature_df


def IMR90(feature_cols):
	#traing_df = pd.read_csv('../data/IMR90/midfile/IMR90_pairs.csv')
	#stable_required_df = traing_df
	stable_required_df = pd.read_csv('../data/IMR90/muti-midfile/train.csv')
	stable_required_df = stable_required_df[use_cols]
	
	all_chrome_seq_feature = pd.DataFrame()

	chr_names = ['chr1','chr2','chr3','chr4',
				'chr5','chr6','chr7','chr8',
				'chr9','chr10','chr11','chr12',
				'chr13','chr14','chr15','chr16',
				'chr17','chr18','chr19','chr20',
				'chr21','chr22','chrX','chrY']
	
	for chr_name in chr_names:	
		required_info_df = stable_required_df[stable_required_df['chrom-Enh']==chr_name]
		if required_info_df.shape[0] == 0:
			continue
		currrent_sequence_feature = generate_feature(filepath="../data/sequence/hg19.id_" + chr_name + ".fa", 
									feature_cols=feature_cols,
									required_info_df=required_info_df,
									chr_name=chr_name)
		#currrent_sequence_feature.to_csv('./data/IMR90/slices/'+chr_name+'.csv')
		print("process chrome: {} successfull,shape is {}".format(chr_name,currrent_sequence_feature.shape))
		
		all_chrome_seq_feature = all_chrome_seq_feature.append(currrent_sequence_feature)

	#all_chrome_seq_feature = all_chrome_seq_feature.fillna(0)
	#print(all_chrome_seq_feature)
	all_chrome_seq_feature.to_csv('../data/EP-lower/low_seq_feature_imr90.csv',index=False)
	
	print('IMR90 low sequence feature saved')


def K562(feature_cols):
	stable_required_df = pd.read_csv('../data/k562/muti-midfile/train.csv')
	stable_required_df = stable_required_df[use_cols]

	all_chrome_seq_feature = pd.DataFrame()

	chr_names = ['chr1','chr2','chr3','chr4',
				'chr5','chr6','chr7','chr8',
				'chr9','chr10','chr11','chr12',
				'chr13','chr14','chr15','chr16',
				'chr17','chr18','chr19','chr20',
				'chr21','chr22','chrX','chrY']
	
	for chr_name in chr_names:	
		required_info_df = stable_required_df[stable_required_df['chrom-Enh']==chr_name]
		if required_info_df.shape[0] == 0:
			continue
		currrent_sequence_feature = generate_feature(filepath="../data/sequence/hg19.id_" + chr_name + ".fa", 
									feature_cols=feature_cols,
									required_info_df=required_info_df,
									chr_name=chr_name)
		#currrent_sequence_feature.to_csv('./data/IMR90/slices/'+chr_name+'.csv')
		print("process chrome: {} successfull,shape is {}".format(chr_name,currrent_sequence_feature.shape))
		
		all_chrome_seq_feature = all_chrome_seq_feature.append(currrent_sequence_feature)

	#all_chrome_seq_feature = all_chrome_seq_feature.fillna(0)
	#print(all_chrome_seq_feature)
	all_chrome_seq_feature.to_csv('../data/EP-lower/low_seq_feature_K562.csv',index=False)
	
	print('k562 low sequence feature saved')



def GM12878(feature_cols):
	stable_required_df = pd.read_csv('../data/GM12878/muti-midfile/train.csv')
	stable_required_df = stable_required_df[use_cols]
	
	all_chrome_seq_feature = pd.DataFrame()

	chr_names = ['chr1','chr2','chr3','chr4',
				'chr5','chr6','chr7','chr8',
				'chr9','chr10','chr11','chr12',
				'chr13','chr14','chr15','chr16',
				'chr17','chr18','chr19','chr20',
				'chr21','chr22','chrX','chrY']
	
	for chr_name in chr_names:	
		required_info_df = stable_required_df[stable_required_df['chrom-Enh']==chr_name]
		if required_info_df.shape[0] == 0:
			continue
		currrent_sequence_feature = generate_feature(filepath="../data/sequence/hg19.id_" + chr_name + ".fa", 
									feature_cols=feature_cols,
									required_info_df=required_info_df,
									chr_name=chr_name)
		#currrent_sequence_feature.to_csv('./data/IMR90/slices/'+chr_name+'.csv')
		print("process chrome: {} successfull,shape is {}".format(chr_name,currrent_sequence_feature.shape))
		
		all_chrome_seq_feature = all_chrome_seq_feature.append(currrent_sequence_feature)

	#all_chrome_seq_feature = all_chrome_seq_feature.fillna(0)
	#print(all_chrome_seq_feature)
	all_chrome_seq_feature.to_csv('../data/EP-lower/low_seq_feature_GM12878.csv',index=False)
	
	print('GM12878 low sequence feature saved')

def HelaS3(feature_cols):
	stable_required_df = pd.read_csv('../data/Hela-S3/muti-midfile/train.csv')
	stable_required_df = stable_required_df[use_cols]

	all_chrome_seq_feature = pd.DataFrame()

	chr_names = ['chr1','chr2','chr3','chr4',
				'chr5','chr6','chr7','chr8',
				'chr9','chr10','chr11','chr12',
				'chr13','chr14','chr15','chr16',
				'chr17','chr18','chr19','chr20',
				'chr21','chr22','chrX','chrY']
	
	for chr_name in chr_names:	
		required_info_df = stable_required_df[stable_required_df['chrom-Enh']==chr_name]
		if required_info_df.shape[0] == 0:
			continue
		currrent_sequence_feature = generate_feature(filepath="../data/sequence/hg19.id_" + chr_name + ".fa", 
									feature_cols=feature_cols,
									required_info_df=required_info_df,
									chr_name=chr_name)
		#currrent_sequence_feature.to_csv('./data/IMR90/slices/'+chr_name+'.csv')
		print("process chrome: {} successfull,shape is {}".format(chr_name,currrent_sequence_feature.shape))
		
		all_chrome_seq_feature = all_chrome_seq_feature.append(currrent_sequence_feature)

	#all_chrome_seq_feature = all_chrome_seq_feature.fillna(0)
	#print(all_chrome_seq_feature)
	all_chrome_seq_feature.to_csv('../data/EP-lower/low_seq_feature_Hela-S3.csv',index=False)
	
	print('Hela-S3 low sequence feature saved')


def HUVEC(feature_cols):
	stable_required_df = pd.read_csv('../data/HUVEC/muti-midfile/train.csv')
	stable_required_df = stable_required_df[use_cols]
	
	all_chrome_seq_feature = pd.DataFrame()

	chr_names = ['chr1','chr2','chr3','chr4',
				'chr5','chr6','chr7','chr8',
				'chr9','chr10','chr11','chr12',
				'chr13','chr14','chr15','chr16',
				'chr17','chr18','chr19','chr20',
				'chr21','chr22','chrX','chrY']
	
	for chr_name in chr_names:	
		required_info_df = stable_required_df[stable_required_df['chrom-Enh']==chr_name]
		if required_info_df.shape[0] == 0:
			continue
		currrent_sequence_feature = generate_feature(filepath="../data/sequence/hg19.id_" + chr_name + ".fa", 
									feature_cols=feature_cols,
									required_info_df=required_info_df,
									chr_name=chr_name)
		#currrent_sequence_feature.to_csv('./data/IMR90/slices/'+chr_name+'.csv')
		print("process chrome: {} successfull,shape is {}".format(chr_name,currrent_sequence_feature.shape))
		
		all_chrome_seq_feature = all_chrome_seq_feature.append(currrent_sequence_feature)

	#all_chrome_seq_feature = all_chrome_seq_feature.fillna(0)
	#print(all_chrome_seq_feature)
	all_chrome_seq_feature.to_csv('../data/EP-lower/low_seq_feature_HUVEC.csv',index=False)
	
	print('HUVEC low sequence feature saved')


def NHEK(feature_cols):
	stable_required_df = pd.read_csv('../data/NHEK/muti-midfile/train.csv')


	stable_required_df = stable_required_df[use_cols]

	all_chrome_seq_feature = pd.DataFrame()

	chr_names = ['chr1','chr2','chr3','chr4',
				'chr5','chr6','chr7','chr8',
				'chr9','chr10','chr11','chr12',
				'chr13','chr14','chr15','chr16',
				'chr17','chr18','chr19','chr20',
				'chr21','chr22','chrX','chrY']
	
	for chr_name in chr_names:	
		required_info_df = stable_required_df[stable_required_df['chrom-Enh']==chr_name]
		if required_info_df.shape[0] == 0:
			continue
		currrent_sequence_feature = generate_feature(filepath="../data/sequence/hg19.id_" + chr_name + ".fa", 
									feature_cols=feature_cols,
									required_info_df=required_info_df,
									chr_name=chr_name)
		#currrent_sequence_feature.to_csv('./data/IMR90/slices/'+chr_name+'.csv')
		print("process chrome: {} successfull,shape is {}".format(chr_name,currrent_sequence_feature.shape))
		
		all_chrome_seq_feature = all_chrome_seq_feature.append(currrent_sequence_feature,ignore_index=True)

	#all_chrome_seq_feature = all_chrome_seq_feature.fillna(0)
	#print(all_chrome_seq_feature)
	all_chrome_seq_feature.to_csv('../data/EP-lower/low_seq_feature_NHEK.csv',index=False)
	
	print('NHEK low sequence feature saved')



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='extract varies cellline features for is it in TAD')
	parser.add_argument("--cellline",type=str,help="select celle line, IMR90 GM12878 K562")
	parser.add_argument("--cols",type=str,help="select feature extract region,window enhancer full promoter2k")
	args = parser.parse_args()

	cellline = args.cellline
	cols = args.cols

	if cellline == 'K562':
		K562(int(cols))
	elif cellline =='GM12878':
		GM12878(int(cols))
	elif cellline == 'IMR90':
		IMR90(int(cols))
	elif cellline == 'HUVEC':
		HUVEC(int(cols))
	elif cellline == 'Hela-S3':
		HelaS3(int(cols))
	elif cellline == 'NHEK':
		NHEK(int(cols))
	else:
		print('wrong argument')




