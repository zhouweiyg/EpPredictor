import pandas as pd
import numpy as np
import random
import os

from sklearn.utils import shuffle

EP_COLUMNS_NAME = ['chrom-Enh','chromStart','chromEnd','Gene ID','chrom-Gene','TSS','Transcript',
					'signalValue','EP Score']


TAD_COLUMNS_NAME = ['chrom','tad_start','tad_end']

BAD_COLUMNS_NAME = ['chrom','chromStart','chromEnd','name',
					'score','strand','signalValue','pValue',
					'qValue','peak']

def load_tad(tad_fpath):
	tad_df = pd.read_csv(tad_fpath, sep='\t', names=TAD_COLUMNS_NAME)
	return tad_df

def load_ep(ep_fpath,sep):
	ep_df = pd.read_csv(ep_fpath, sep=sep, names=EP_COLUMNS_NAME)
	return ep_df

def load_chrom():
	chrom_vec = ['chr'+str(x) for x in range(1,25)]
	chrom_vec.extend(['chrX','chrY'])
	return chrom_vec

def generate_enhancer_promoter_pair(ep_df):
	"""

	"""

	std_ep_pair = ep_df[['chrom-Enh','chromStart','chromEnd','TSS']]
	min_ep_gap = abs((std_ep_pair['chromEnd']-std_ep_pair['chromStart']).min())
	max_ep_gap = abs((std_ep_pair['chromEnd']-std_ep_pair['chromStart']).max())

	fake_samples = []
	for enhancer in std_ep_pair[['chrom-Enh','chromStart','chromEnd']].values:
		for promoter in std_ep_pair['TSS'].values:
			gap = abs(enhancer[-1]-promoter)
			if gap>min_ep_gap and gap<max_ep_gap:
				current_sample = np.r_[enhancer, promoter]
				fake_samples.append(current_sample)
	fake_samples = random.sample(fake_samples, std_ep_pair.shape[0])
	fake_ep_pair = pd.DataFrame(fake_samples, columns=['chrom-Enh','chromStart','chromEnd','TSS'])

	return std_ep_pair, fake_ep_pair


def load_pairs(hic_ep_file):
	use_cols = ['chrom-Enh','chromStart','chromEnd','TSS', 'label']
	samples = pd.read_csv(hic_ep_file)#'../midfile/hic_ep.csv'
	samples = samples[use_cols]
	'''这部份是没有加hic时候的标签样本
	real_samples = pd.read_csv('../midfile/std_ep_pair.csv')
	fake_samples = pd.read_csv('../midfile/fake_ep_pair.csv')

	real_samples['label'] = 1
	fake_samples['label'] = 0

	samples = pd.concat([real_samples, fake_samples], ignore_index=True)
	'''
	return samples


def load_targetfinder_pairs(pairs_path):
	use_cols = ['enhancer_chrom','enhancer_start','enhancer_end','promoter_start','label']
	samples = pd.read_csv(pairs_path)
	samples = samples[use_cols]

	samples = samples.rename(columns={'enhancer_chrom':'chrom-Enh',
									'enhancer_start':'chromStart',
									'enhancer_end':'chromEnd',
									'promoter_start':'TSS',
									})

	samples_1_nums = samples[samples['label']==1].shape[0]
	samples_0_nums = samples[samples['label']==0].shape[0]

	min_samples_caterogy_num = min(samples_1_nums,samples_0_nums)
	samples_1 = shuffle(samples[samples['label']==1])
	samples_0 = shuffle(samples[samples['label']==0])

	samples_1 = samples_1[:min_samples_caterogy_num]
	samples_0 = samples_0[:min_samples_caterogy_num]
	samples = samples_1.append(samples_0,ignore_index=True)

	return samples


def check_feature_extracted(histone_name,transcript_root_path):
	transcript_root_path = transcript_root_path.split('/')
	
	feature_path = '/home/special/user/dynamic/zhouyc/ProjectBiotechnology/EP-interaction/data/'+transcript_root_path[2]+'/muti-midfile/enhancerfeatures/'
	extracted_features = os.listdir(feature_path)
	
	for feature_file in extracted_features:
		feature_file = feature_file.split('_')[0]
		if feature_file == histone_name:
			return True
	return False


def load_histons(transcript_root_path, histone_modification_root_path):
	"""
	histone_modification_root_path: fail
	"""
	histone_dict = {}

	bind_histons_path = os.listdir(transcript_root_path)#'../bindhistone/'
	for bindhiston_name in bind_histons_path:
		abs_path = transcript_root_path+bindhiston_name
		name = bindhiston_name.strip('.bed')
		if(check_feature_extracted(name,transcript_root_path)):
			print("{} has been extracted, continue".format(name))
			continue
		bindhiston_df = pd.read_csv(abs_path, names=BAD_COLUMNS_NAME, sep='\t')
		histone_dict[name] = bindhiston_df
	return histone_dict
