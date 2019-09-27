
import pandas as pd
import os
import numpy as np


def merge(enhancer_features_path, window_features_path, 
		  patch_features_path,promoter2k_features_path,
		  complete_feature_save_path):
	enhancer_features_path_list = os.listdir(enhancer_features_path)
	window_features_path_list = os.listdir(window_features_path)
	patch_features_path_list = os.listdir(patch_features_path)
	promoter2k_features_path_list = os.listdir(promoter2k_features_path)

	abs_paths = []
	for enhancer_path in enhancer_features_path_list:
		abs_enhancer_path = enhancer_features_path+enhancer_path
		abs_paths.append(abs_enhancer_path)

	for window_path in window_features_path_list:
		abs_window_path = window_features_path+window_path
		abs_paths.append(abs_window_path)

	for patch_path in patch_features_path_list:
		abs_patch_path = patch_features_path+patch_path
		abs_paths.append(abs_patch_path)

	for promoter2k_path in promoter2k_features_path_list:
		abs_promoter2k_path =  promoter2k_features_path+promoter2k_path
		abs_paths.append(abs_promoter2k_path)

	complete_data = pd.read_csv(abs_paths[0])
	print(complete_data.shape)
	for abs_path in abs_paths[1:]:
		current_df = pd.read_csv(abs_path)
		#current_df = current_df.drop(['chrom-Enh','chromStart','chromEnd','TSS','label','id'],axis=1)
		#'''
		complete_data = complete_data.merge(current_df, how='left', on=['chrom-Enh',
																		'chromStart',
																		'chromEnd',
																		'TSS',
																		'label',
																		'id'])
		#'''
		#complete_data = pd.concat([complete_data,current_df],axis=1)
		print(complete_data.shape,abs_path)
	
	complete_data.to_csv(complete_feature_save_path, index=False)
	print(complete_data.shape)

def merge_enhancer(use_real_ep_folder, use_fake_ep_folder, 
				  cellline, region):
	use_real_ep_files = os.listdir(use_real_ep_folder)
	use_fake_ep_files = os.listdir(use_fake_ep_folder)
	intersect = np.intersect1d(use_real_ep_files, use_fake_ep_files)

	all_merged = pd.DataFrame()
	for file in intersect:
		histon_name = file.split('_')[0]
		if region == 'enhancerfeatures':		
			use_real_ep_columns = ['chrom-Enh','chromStart','chromEnd','TSS','label','id']
			enhancer_features_cols = [histon_name+'_'+x for x in ['enhancer_avg_sigval','enhancer_max_sigval','enhancer_min_sigval',
															 'enhancer_std_sigval','enhancer_var_sigval','enhancer_avg_peak',
															 'enhancer_max_peak','enhancer_min_peak','enhancer_std_peak','enhancer_var_peak',
													 		  'enhancer_peak_density']]	
			use_real_ep_columns.extend(enhancer_features_cols)

			real_data = pd.read_csv(use_real_ep_folder+file,names=use_real_ep_columns)
		else:
			real_data = pd.read_csv(use_real_ep_folder+file)
		fake_data = pd.read_csv(use_fake_ep_folder+file)
		merged_tmp = fake_data.append(real_data)
		merged_tmp.to_csv('../merged-ef/'+cellline+'/'+region+'/'+file,index=False)
		print(merged_tmp.shape)

if __name__ == '__main__':
	'''
	merge_enhancer('../../data/GM12878/midfile/enhancerfeatures/','../data/GM12878/muti-midfile/enhancerfeatures/',
				'GM12878','enhancerfeatures')


	merge_enhancer('../../data/GM12878/midfile/windowfeatures/','../data/GM12878/muti-midfile/windowfeatures/',
				'GM12878','windowfeatures')
	merge_enhancer('../../data/GM12878/midfile/patchfeatures/','../data/GM12878/muti-midfile/patchfeatures/',
				'GM12878','patchfeatures')	


	merge_enhancer('../../data/GM12878/midfile/promoter2k/','../data/GM12878/muti-midfile/promoter2k/',
				'GM12878','promoter2k')	

	'''

	'''
	merge('../data/NHEK/muti-midfile/enhancerfeatures/',
		'../data/NHEK/muti-midfile/windowfeatures/',
		'../data/NHEK/muti-midfile/patchfeatures/',
		'../data/NHEK/muti-midfile/promoter2k/',
		'../data/NHEK/muti-midfile/train.csv')
	'''
	for item in ['NHEK']:#['k562','GM12878','IMR90','HUVEC','Hela-S3','NHEK']:
		merge('../data/{}/muti-midfile/enhancerfeatures/'.format(item),
			'../data/{}/muti-midfile/windowfeatures/'.format(item),
			'../data/{}/muti-midfile/patchfeatures/'.format(item),
			'../data/{}/muti-midfile/promoter2k/'.format(item),
			'../data/{}/muti-midfile/train.csv'.format(item))
