import pandas as pd
import numpy as np
import pandas as pd
from . import utils


import threading

class ThreadNormalEnhancerFeature(threading.Thread):
	def __init__(self, chrom_vec, samples, histon_df, histon_name, feature_save_path):
		"""
		@samples pandas like [['id,chrom-Enh,chromStart,chromEnd,TSS,label']]
		"""
		threading.Thread.__init__(self)
		self.chrom_vec = chrom_vec
		self.samples = samples
		self.histone_df = histon_df
		self.histon_name = histon_name
		self.feature_save_path = feature_save_path

	def extract_promoter(self):
		pass


	def run(self):
		enhancer_features = self.samples
		histon_df =self.histone_df
		original_cols =  self.samples.columns.tolist()

		features_cols = [self.histon_name+'_'+x for x in ['avg_sigval','max_sigval','min_sigval','std_sigval','var_sigval',
													 'avg_peak','max_peak','min_peak','std_peak','var_peak',
													 'peak_density']]

		features_cols.insert(0,'id')
		features_diff_histone = []
	
		for sample in self.samples.values:
			chrom_name = sample[0]
			chrom_start = sample[1]
			chrom_end = sample[2]
			idx = sample[-1]
				
			subset = histon_df[(histon_df['chrom']==chrom_name) &\
									(histon_df['chromStart']>=chrom_start) &\
									(histon_df['chromEnd']<=chrom_end)]
				
			avg_sigval= subset['signalValue'].mean()
			max_sigval = subset['signalValue'].max()
			min_sigval = subset['signalValue'].min()
			std_sigval = subset['signalValue'].std()
			var_sigval = subset['signalValue'].var()

			avg_peak = subset['peak'].mean()
			max_peak = subset['peak'].max()
			min_peak = subset['peak'].min()
			std_peak = subset['peak'].std()
			var_peak = subset['peak'].var()

			peak_density = subset.shape[0]

			
			


			features = [idx, avg_sigval, max_sigval, min_sigval, std_sigval, var_sigval,
						avg_peak, max_peak, min_peak, std_peak, var_peak,
						peak_density]

	
			features_diff_histone.append(features)
		
		features_diff_histone = pd.DataFrame(features_diff_histone, columns=features_cols)
		enhancer_features = enhancer_features.merge(features_diff_histone, on='id',how='left')
		
		enhancer_features.to_csv(self.feature_save_path+self.histon_name+'_enhancer_features.csv',index=False)
		print(self.histon_name, enhancer_features.shape)




class ThreadNormalWondowrFeature(threading.Thread):
	def __init__(self, chrom_vec, samples, histon_df, histon_name, feature_save_path):
		"""
		@samples pandas like [['id,chrom-Enh,chromStart,chromEnd,TSS,label']]
		"""
		threading.Thread.__init__(self)
		self.chrom_vec = chrom_vec
		self.samples = samples
		self.histone_df = histon_df
		self.histon_name = histon_name
		self.feature_save_path = feature_save_path


	def run(self):
		window_features = self.samples
		histon_df =self.histone_df
		original_cols =  self.samples.columns.tolist()

		features_cols = [self.histon_name+'_'+x for x in ['win_avg_sigval','win_max_sigval','win_min_sigval','win_std_sigval','win_var_sigval',
													 'win_avg_peak','win_max_peak','win_min_peak','win_std_peak','win_var_peak',
													 'win_peak_density']]

		features_cols.insert(0,'id')
		features_diff_histone = []
	
		for sample in self.samples.values:
			chrom_name = sample[0]
			chrom_start = sample[1]
			chrom_end = sample[2]
			promoter_loc = sample[3]
			idx = sample[-1]
			

			if promoter_loc>chrom_end:
				subset = histon_df[(histon_df['chrom']==chrom_name) &\
										(histon_df['chromStart']>=chrom_end) &\
										(histon_df['chromEnd']<=promoter_loc)]
			else:
				subset = histon_df[(histon_df['chrom']==chrom_name) &\
									(histon_df['chromStart']>=promoter_loc) &\
									(histon_df['chromEnd']<=chrom_end)]
				
			avg_sigval= subset['signalValue'].mean()
			max_sigval = subset['signalValue'].max()
			min_sigval = subset['signalValue'].min()
			std_sigval = subset['signalValue'].std()
			var_sigval = subset['signalValue'].var()

			avg_peak = subset['peak'].mean()
			max_peak = subset['peak'].max()
			min_peak = subset['peak'].min()
			std_peak = subset['peak'].std()
			var_peak = subset['peak'].var()

			peak_density = subset.shape[0]

			features = [idx, avg_sigval, max_sigval, min_sigval, std_sigval, var_sigval,
						avg_peak, max_peak, min_peak, std_peak, var_peak,
						peak_density]

			features_diff_histone.append(features)
		
		features_diff_histone = pd.DataFrame(features_diff_histone, columns=features_cols)
		window_features = window_features.merge(features_diff_histone, on='id',how='left')
		
		window_features.to_csv(self.feature_save_path+self.histon_name+'_window_features.csv',index=False)
		print(self.histon_name, window_features.shape)


class ThreadPatchFeature(threading.Thread):
	def __init__(self, chrom_vec, samples, histon_df, histon_name, feature_save_path):
		"""
		@samples pandas like [['id,chrom-Enh,chromStart,chromEnd,TSS,label']]
		"""
		threading.Thread.__init__(self)
		self.chrom_vec = chrom_vec
		self.samples = samples
		self.histone_df = histon_df
		self.histon_name = histon_name
		self.feature_save_path = feature_save_path


	def run(self):
		patch_features = self.samples
		histon_df =self.histone_df
		original_cols =  self.samples.columns.tolist()

		features_cols = [self.histon_name+'_'+x for x in ['patch_avg_sigval','patch_max_sigval','patch_min_sigval','patch_std_sigval','patch_var_sigval',
													 'patch_avg_peak','patch_max_peak','patch_min_peak','patch_std_peak','patch_var_peak',
													 'patch_peak_density']]

		features_cols.insert(0,'id')
		features_diff_histone = []
	
		for sample in self.samples.values:
			chrom_name = sample[0]
			chrom_start = sample[1]
			chrom_end = sample[2]
			promoter_loc = sample[3]
			idx = sample[-1]
			

			if promoter_loc>chrom_end:
				subset = histon_df[(histon_df['chrom']==chrom_name) &\
										(histon_df['chromStart']>=chrom_start) &\
										(histon_df['chromEnd']<=promoter_loc)]
			else:
				subset = histon_df[(histon_df['chrom']==chrom_name) &\
									(histon_df['chromStart']>=chrom_start) &\
									(histon_df['chromEnd']<=chrom_end)]
				
			avg_sigval= subset['signalValue'].mean()
			max_sigval = subset['signalValue'].max()
			min_sigval = subset['signalValue'].min()
			std_sigval = subset['signalValue'].std()
			var_sigval = subset['signalValue'].var()

			avg_peak = subset['peak'].mean()
			max_peak = subset['peak'].max()
			min_peak = subset['peak'].min()
			std_peak = subset['peak'].std()
			var_peak = subset['peak'].var()

			peak_density = subset.shape[0]

			features = [idx, avg_sigval, max_sigval, min_sigval, std_sigval, var_sigval,
						avg_peak, max_peak, min_peak, std_peak, var_peak,
						peak_density]
			
			features_diff_histone.append(features)
		
		features_diff_histone = pd.DataFrame(features_diff_histone, columns=features_cols)
		patch_features = patch_features.merge(features_diff_histone, on='id',how='left')
		
		patch_features.to_csv(self.feature_save_path+self.histon_name+'_patch_features.csv',index=False)
		print(self.histon_name, patch_features.shape)


class ThreadPromoterDistanceFeature(threading.Thread):
	def __init__(self, chrom_vec, samples, histon_df, histon_name, feature_save_path):
		"""
		@samples pandas like [['id,chrom-Enh,chromStart,chromEnd,TSS,label']]
		"""
		threading.Thread.__init__(self)
		self.chrom_vec = chrom_vec
		self.samples = samples
		self.histone_df = histon_df
		self.histon_name = histon_name
		self.feature_save_path = feature_save_path


	def run(self):
		promoter2k_features = self.samples
		histon_df =self.histone_df
		original_cols =  self.samples.columns.tolist()

		features_cols = [self.histon_name+'_'+x for x in ['promoter_2k_avg_sigval','promoter_2k_max_sigval',
														   'promoter_2k_min_sigval','promoter_2k_std_sigval',
														   'promoter_2k_var_sigval','promoter_2k_avg_peak',
														   'promoter_2k_max_peak','promoter_2k_min_peak',
														   'promoter_2k_std_peak','promoter_2k_var_peak',
													 		'promoter_2k_peak_density']]

		features_cols.insert(0,'id')
		features_diff_histone = []
	
		for sample in self.samples.values:
			chrom_name = sample[0]
			chrom_start = sample[1]
			chrom_end = sample[2]
			promoter_loc = sample[3]

			positive_promoter_loc = promoter_loc+2000
			negative_promoter_loc = promoter_loc-2000
			idx = sample[-1]
			

			if promoter_loc>chrom_end:
				subset = histon_df[(histon_df['chrom']==chrom_name)&\
								(histon_df['chromStart'])>=negative_promoter_loc&\
								(histon_df['chromEnd']<=positive_promoter_loc)]
			
			elif promoter_loc>chrom_start and promoter_loc<chrom_end:
				negative_gap = promoter_loc - chrom_start-2000
				positive_gap = chrom_end - promoter_loc-2000
				subset = histon_df[(histon_df['chrom']==chrom_name) &\
									((histon_df['chromStart']+negative_gap)>=chrom_start) &\
									((histon_df['chromEnd']-positive_gap)<=chrom_end)]
			else:
				subset = histon_df[(histon_df['chrom']==chrom_name)&\
								(histon_df['chromStart'])>=negative_promoter_loc&\
								(histon_df['chromEnd']<=positive_promoter_loc)]
				


			
			avg_sigval= subset['signalValue'].mean()
			max_sigval = subset['signalValue'].max()
			min_sigval = subset['signalValue'].min()
			std_sigval = subset['signalValue'].std()
			var_sigval = subset['signalValue'].var()

			avg_peak = subset['peak'].mean()
			max_peak = subset['peak'].max()
			min_peak = subset['peak'].min()
			std_peak = subset['peak'].std()
			var_peak = subset['peak'].var()

			peak_density = subset.shape[0]

			features = [idx, avg_sigval, max_sigval, min_sigval, std_sigval, var_sigval,
						avg_peak, max_peak, min_peak, std_peak, var_peak,
						peak_density]
			
			features_diff_histone.append(features)
		
		features_diff_histone = pd.DataFrame(features_diff_histone, columns=features_cols)
		promoter2k_features = promoter2k_features.merge(features_diff_histone, on='id',how='left')
		
		promoter2k_features.to_csv(self.feature_save_path+self.histon_name+'_promoter2k_features.csv',index=False)
		print(self.histon_name, promoter2k_features.shape)
