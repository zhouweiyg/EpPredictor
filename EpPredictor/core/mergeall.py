
import pandas as pd
import os


'''
enhancer_features_path = '../midfile/enhancerfeatures/'
window_features_path = '../midfile/windowfeatures/'
patch_features_path = '../midfile/patchfeatures/'
'''

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
		complete_data = complete_data.merge(current_df, how='left', on=['chrom-Enh',
																		'chromStart',
																		'chromEnd',
																		'TSS',
																		'label',
																		'id'])
		print(complete_data.shape)
	print(complete_data.shape)
	complete_data.to_csv(complete_feature_save_path, index=False)