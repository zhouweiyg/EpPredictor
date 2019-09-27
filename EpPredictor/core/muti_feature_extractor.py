from multiprocessing import cpu_count
import utils
import pandas as pd
from multiprocessing import Process
import argparse

CPU_CORE_NUMS = cpu_count()


def subset_featuras(subset,idx):
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

	return features

def feature_processor(chrom_vec, samples, histon_df_dict, 
								histon_name_vec,
								feature_save_path):
	
	for histon_name in histon_name_vec:
		enhancer_features = samples
		original_cols =  samples.columns.tolist()

		histon_df = histon_df_dict[histon_name]
		enhancer_features_cols = [histon_name+'_'+x for x in ['enhancer_avg_sigval','enhancer_max_sigval','enhancer_min_sigval',
															 'enhancer_std_sigval','enhancer_var_sigval','enhancer_avg_peak',
															 'enhancer_max_peak','enhancer_min_peak','enhancer_std_peak','enhancer_var_peak',
													 		  'enhancer_peak_density']]

		window_features_cols = [histon_name+'_'+x for x in ['win_avg_sigval','win_max_sigval','win_min_sigval','win_std_sigval','win_var_sigval',
													 'win_avg_peak','win_max_peak','win_min_peak','win_std_peak','win_var_peak',
													 'win_peak_density']]

		patch_features_cols = [histon_name+'_'+x for x in ['patch_avg_sigval','patch_max_sigval','patch_min_sigval','patch_std_sigval','patch_var_sigval',
													 'patch_avg_peak','patch_max_peak','patch_min_peak','patch_std_peak','patch_var_peak',
													 'patch_peak_density']]

		promoter2k_features_cols = [histon_name+'_'+x for x in ['promoter_2k_avg_sigval','promoter_2k_max_sigval',
														   'promoter_2k_min_sigval','promoter_2k_std_sigval',
														   'promoter_2k_var_sigval','promoter_2k_avg_peak',
														   'promoter_2k_max_peak','promoter_2k_min_peak',
														   'promoter_2k_std_peak','promoter_2k_var_peak',
													 		'promoter_2k_peak_density']]

		enhancer_features_cols.insert(0,'id')
		window_features_cols.insert(0,'id')
		patch_features_cols.insert(0,'id')
		promoter2k_features_cols.insert(0,'id')

		enhancer_features_diff_histone = []
		window_features_diff_histone = []
		promoter2k_features_diff_histone = []
		patch_features_diff_histone = []
		
		for sample in samples.values:
			chrom_name = sample[0]
			chrom_start = sample[1]
			chrom_end = sample[2]
			promoter_loc = sample[3]
			idx = sample[-1]
			
			positive_promoter_loc = promoter_loc+2000
			negative_promoter_loc = promoter_loc-2000

			##enhancer
			subset_enhancer = histon_df[(histon_df['chrom']==chrom_name) &\
								(histon_df['chromStart']>=chrom_start) &\
								(histon_df['chromEnd']<=chrom_end)]


			#window
			if promoter_loc>chrom_end:
				subset_window = histon_df[(histon_df['chrom']==chrom_name) &\
										(histon_df['chromStart']>=chrom_end) &\
										(histon_df['chromEnd']<=promoter_loc)]
			else:
				subset_window = histon_df[(histon_df['chrom']==chrom_name) &\
									(histon_df['chromStart']>=promoter_loc) &\
									(histon_df['chromEnd']<=chrom_end)]
			
			#########promoter2k
			if promoter_loc>chrom_end:
				subset_promotet2k = histon_df[(histon_df['chrom']==chrom_name)&\
								(histon_df['chromStart'])>=negative_promoter_loc&\
								(histon_df['chromEnd']<=positive_promoter_loc)]
			
			elif promoter_loc>chrom_start and promoter_loc<chrom_end:
				negative_gap = promoter_loc - chrom_start-2000
				positive_gap = chrom_end - promoter_loc-2000
				subset_promotet2k = histon_df[(histon_df['chrom']==chrom_name) &\
									((histon_df['chromStart']+negative_gap)>=chrom_start) &\
									((histon_df['chromEnd']-positive_gap)<=chrom_end)]
			else:
				subset_promotet2k = histon_df[(histon_df['chrom']==chrom_name)&\
								(histon_df['chromStart'])>=negative_promoter_loc&\
								(histon_df['chromEnd']<=positive_promoter_loc)]
			
			####patch
			if promoter_loc>chrom_end:
				subset_patch = histon_df[(histon_df['chrom']==chrom_name) &\
										(histon_df['chromStart']>=chrom_start) &\
										(histon_df['chromEnd']<=promoter_loc)]
			else:
				subset_patch = histon_df[(histon_df['chrom']==chrom_name) &\
									(histon_df['chromStart']>=chrom_start) &\
									(histon_df['chromEnd']<=chrom_end)]
			
			sub_enhancer_features = subset_featuras(subset_enhancer,idx)
			sub_window_features = subset_featuras(subset_window,idx)
			sub_patch_features = subset_featuras(subset_window,idx)
			sub_promoter2k_features = subset_featuras(subset_promotet2k,idx)
			
			enhancer_features_diff_histone.append(sub_enhancer_features)
			window_features_diff_histone.append(sub_window_features)
			patch_features_diff_histone.append(sub_patch_features)
			promoter2k_features_diff_histone.append(sub_promoter2k_features)
			
		enhancer_features_diff_histone = pd.DataFrame(enhancer_features_diff_histone, columns=enhancer_features_cols)
		enhancer_features = samples.merge(enhancer_features_diff_histone, on='id',how='left')
		
		window_features_diff_histone = pd.DataFrame(window_features_diff_histone, columns=window_features_cols)
		window_features = samples.merge(window_features_diff_histone, on='id',how='left')
		
		patch_features_diff_histone = pd.DataFrame(patch_features_diff_histone, columns=patch_features_cols)
		patch_features = samples.merge(patch_features_diff_histone, on='id',how='left')
		
		promoter2k_features_diff_histone = pd.DataFrame(promoter2k_features_diff_histone, columns=promoter2k_features_cols)
		promoter2k_features = samples.merge(promoter2k_features_diff_histone, on='id',how='left')
		

		enhancer_features.to_csv(feature_save_path+'enhancerfeatures/'+histon_name+'_enhancer_features.csv',index=False)
		window_features.to_csv(feature_save_path+'windowfeatures/'+histon_name+'_window_features.csv',index=False)
		patch_features.to_csv(feature_save_path+'patchfeatures/'+histon_name+'_patch_features.csv',index=False)
		promoter2k_features.to_csv(feature_save_path+'promoter2k/'+histon_name+'_promoter2k_features.csv',index=False)
		
		print(histon_name, enhancer_features.shape,window_features.shape,patch_features.shape,promoter2k_features.shape)
		

def run_feature_processor():
	processes = list()
	histone_df_dict = utils.load_histons(transcript_root_path='../data/IMR90/transcription factor/', 
										  histone_modification_root_path=None)
	chrom_vec = utils.load_chrom()
	histone_name_vec = list(histone_df_dict.keys())
	#histone_name_vec = histone_name_vec[:35]
	#print(histone_name_vec)
	samples = utils.load_pairs(hic_ep_file='../EP-interaction/data/IMR90/midfile/hic_ep.csv')
	samples = samples[samples['label']==0]
	samples['id'] = samples.index
	print(samples)
	feature_save_path = '../EP-interaction/data/IMR90/muti-midfile/'

	avg_process_histone = len(histone_name_vec)//CPU_CORE_NUMS

	init_process_histone_idx = 0
	for i in range(CPU_CORE_NUMS):
		print ('Process will start, process histon_df nums={}'.format(init_process_histone_idx))
		p = Process(target=feature_processor, args=(chrom_vec,samples,
        													histone_df_dict,
        													histone_name_vec[init_process_histone_idx:init_process_histone_idx+avg_process_histone],
        													feature_save_path,))
		init_process_histone_idx += avg_process_histone
		if i==CPU_CORE_NUMS:
			avg_process_histone = 100
		print ('Process will start, process histon_df nums={}'.format(init_process_histone_idx))
		p.start()
		processes.append(p)
    
	for p in processes:
		p.join()
	print('Process end.')


def run_feature_processor_targetfinder(transcript_root_path,
										interaction_filepath,
										feature_save_path):
	global CPU_CORE_NUMS

	processes = list()
	histone_df_dict = utils.load_histons(transcript_root_path=transcript_root_path, 
										  histone_modification_root_path=None)

	print('histone nums= ',len(histone_df_dict))
	chrom_vec = utils.load_chrom()
	histone_name_vec = list(histone_df_dict.keys())
	#histone_name_vec = histone_name_vec[:35]
	#print(histone_name_vec)
	samples = utils.load_targetfinder_pairs(pairs_path=interaction_filepath)
	#samples = samples[samples['label']==0]
	samples['id'] = samples.index
	print(samples)
	feature_save_path = feature_save_path

	avg_process_histone = len(histone_name_vec)//CPU_CORE_NUMS
	if len(histone_name_vec) < CPU_CORE_NUMS:
		avg_process_histone = 1
		CPU_CORE_NUMS = len(histone_name_vec)

	init_process_histone_idx = 0
	for i in range(CPU_CORE_NUMS):
		print ('Process will start, process histon_df nums={}'.format(init_process_histone_idx))
		p = Process(target=feature_processor, args=(chrom_vec,samples,
        													histone_df_dict,
        													histone_name_vec[init_process_histone_idx:init_process_histone_idx+avg_process_histone],
        													feature_save_path,))
		init_process_histone_idx += avg_process_histone
		if i==CPU_CORE_NUMS:
			avg_process_histone = 100
		print ('Process will start, process histon_df nums={}'.format(init_process_histone_idx))
		p.start()
		processes.append(p)
    
	for p in processes:
		p.join()
	print('Process end.')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='extract varies cellline features')
	parser.add_argument("--chipseq_path",type=str,
						default='../data/IMR90/',
						help="select celle line, IMR90 GM12878 K562")
	
	parser.add_argument("--interaction_path",
						type=str,
						default='../EP-interaction/data/IMR90/midfile/IMR90_pairs.csv',
						help="select feature extract region,window enhancer full promoter2k")

	parser.add_argument("--feature_save_path",
						type=str,
						default='../EP-interaction/data/IMR90/muti-midfile/',
						help="select feature extract region,window enhancer full promoter2k")
	

	args = parser.parse_args()

	protein_path = args.chipseq_path
	interac_path = args.interaction_path
	save_path = args.feature_save_path

	run_feature_processor_targetfinder(transcript_root_path=protein_path,
										interaction_filepath=interac_path,
										feature_save_path=save_path)


