import pandas as pd
import os
import numpy as np

'''
计算各蛋白的重要性
'''
def merge_protine_info(protine_info_path):
	fa = pd.read_excel('./数据信息—秦婉婷180426.xlsx')
	fb = pd.read_excel('./工作簿1GLLrev.xlsx')
	fa = fa.append(fb)

	fa.to_csv(protine_info_path,index=False)

'''
加载蛋白信息，也就是数据信息
'''
def load_protine_info(protine_info_path):
	if not os.path.exists(protine_info_path):
		merge_protine_info(protine_info_path)

	protine_info_df = pd.read_csv(protine_info_path)
	
	return load_protine_info


def choose_best_protine(features_importance_df, protine_info_df, 
						cell_type):
	
	protine_info_df = protine_info_df[protine_info_df['Biosample ']==cell_type]
	features_importance_df = features_importance_df.sort_values(by='score', ascending=False)
	protine_importance_info = []

	for each_feature in features_importance_df.values:
		protine_name = each_feature[1].split('_')[0]
		feature_score = each_feature[0]
		corr_protine_info = protine_info_df[protine_info_df['Accession']==protine_name]

		Experiment = np.nan; Assay=np.nan; Target=np.nan; Biosample=cell_type; Accession=protine_name
		
		if corr_protine_info.shape[0] == 0:
			tmp_info = [Experiment, Assay, Target, Biosample, Accession, feature_score]
			#protine_importance_info.append(tmp_info)
		else:
			corr_protine_info = corr_protine_info.values[0,:]
			Experiment = corr_protine_info[1]
			Assay=corr_protine_info[2]
			Target=corr_protine_info[3]
			Biosample=cell_type
			Accession=protine_name
			tmp_info = [Experiment, Assay, Target, Biosample, Accession, feature_score]
			protine_importance_info.append(tmp_info)

	protine_importance_info = pd.DataFrame(protine_importance_info, 
											columns=['Experiment', 
											'		Assay', 'Target', 
													'Biosample', 'Accession', 
													'feature_score'])

	protine_importance_info = protine_importance_info.drop_duplicates(subset=['Accession'])
	
	return protine_importance_info

def choose_best_protine2(features_importance_df, protine_info_df, 
						cell_type):
	
	protine_info_df = protine_info_df[protine_info_df['Biosample ']==cell_type]
	features_importance_df = features_importance_df.sort_values(by='score', ascending=False)
	protine_importance_info = []

	for each_feature in features_importance_df.values:
		split_features_name = each_feature[1].split('_')
		tail = "";feature_region = ""
		if split_features_name[1] == 'promoter':
			feature_region = 'promoter2k'
		elif split_features_name[1] == 'win':
			feature_region = 'window'
		elif split_features_name[1] == 'patch':
			feature_region = 'full region'
		else:
			feature_region = 'enhancer'

		protine_name = each_feature[1].split('_')[0]
		feature_score = each_feature[0]
		corr_protine_info = protine_info_df[protine_info_df['Accession']==protine_name]
		Experiment = np.nan; Assay=np.nan; Target=np.nan; Biosample=cell_type; Accession=protine_name
		if corr_protine_info.shape[0] == 0:
			tmp_info = [Target,feature_region,feature_score]
			protine_importance_info.append(tmp_info)
		else:
			corr_protine_info = corr_protine_info.values[0,:]
			Experiment = corr_protine_info[1]
			Assay=corr_protine_info[2]
			Target=corr_protine_info[3]+tail
			Biosample=cell_type
			Accession=protine_name
			tmp_info = [Target,feature_region,feature_score]
			protine_importance_info.append(tmp_info)

	protine_importance_info = pd.DataFrame(protine_importance_info, 
											columns=['Target','feature_region','feature_score'])

	protine_importance_info = protine_importance_info.drop_duplicates(subset=['Target','feature_region'])
	
	return protine_importance_info



'''
if __name__ == '__main__':
	features_importance_df = pd.read_csv('../data/k562/midfile/feature_importance_k562_categoty.csv')
	protine_info_df = pd.read_excel('./protine_info.xlsx')
	cell_type = 'Homo sapiens K562'
	cell_protine_imp = choose_best_protine(features_importance_df,protine_info_df, cell_type)
	print(cell_protine_imp)
	cell_protine_imp.to_csv('../data/k562/midfile/k562_imp.csv', index=False)
'''