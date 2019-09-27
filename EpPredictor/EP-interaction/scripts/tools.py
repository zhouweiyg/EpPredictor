
import pandas as pd
import sys
import os
import numpy as np


GTF_FILE_PATH = '../data/gtf/Homo_sapiens.GRCh37.87.chr.gtf'
K562_EP_PATH = '../../data/k562/midfile/hic_ep.csv'
IMR90_EP_PATH = '../../data/IMR90/midfile/hic_ep.csv'
GM12878_EP_PATH = '../../data/GM12878/midfile/hic_ep.csv'



def read_gtf():
	gtf_df = pd.read_csv(GTF_FILE_PATH,sep='\t',skiprows=5,header=None)
	gtf_df = gtf_df.rename(columns={0:'chrom-Enh',1:'method',2:'type',
									3:'start_loc',4:'end_loc',6:'biochain',
									8:'Transcript'})
	rebuild_information = []
	for item in gtf_df['Transcript'].values:
		item = item.split(';')[2]
		item = item.split(' ')[2].strip('"')
		rebuild_information.append(item)
	gtf_df['Transcript'] = rebuild_information

	promoter_loc_vec = []
	for item in gtf_df[['start_loc','end_loc','biochain']].values:
		if item[2] is '+':
			promoter_loc_vec.append(item[0])
		elif item[2] is '-':
			promoter_loc_vec.append(item[1])
		else:
			print("something wrong when judge bio chain")
			sys.exit(1)
	gtf_df['chrom-Enh'] = gtf_df['chrom-Enh'].astype(str)    
	gtf_df['chrom-Enh'] = 'chr'+gtf_df['chrom-Enh']
	gtf_df['TSS'] = promoter_loc_vec
	return gtf_df

def generate_fake_ep(celline):
	
	"""
	根据enhancer-promoter的距离分布
	选择负样本 均值为10000，选择8000~12000距离内的Promoter与已有的enhancer配对
	"""
	def vec2str(vec):
		vec = [str(x) for x in vec]
		vec = " ".join(vec)
		return vec

	gtf_df = read_gtf()
	gtf_df = gtf_df[gtf_df['type']=='transcript']

	save_path = None
	ep_df = None
	if celline == 'K562':
		ep_df = pd.read_csv(K562_EP_PATH,sep=',')
		save_path = '../data/k562/midfile/hic_ep.csv'
	elif celline == 'GM12878':
		ep_df = pd.read_csv(GM12878_EP_PATH, sep=',')
		save_path = '../data/GM12878/midfile/hic_ep.csv'
	elif celline == 'IMR90':
		ep_df = pd.read_csv(IMR90_EP_PATH, sep=',')
		save_path = '../data/IMR90/midfile/hic_ep.csv'
	else:
		print('wrong parameters, parameters should be K562 GM12878 or IMR90')
		sys.exit(1)

	#unique_geneid_vec = np.unique(ep_df['Transcript'])
	transcript_inter_sec = np.intersect1d(gtf_df['Transcript'],ep_df['Transcript'])
	
	sub_gtf_df = gtf_df[~(gtf_df['Transcript'].isin(transcript_inter_sec))]
	sub_gtf_df = sub_gtf_df.sort_values(by=['TSS'],ascending=True).drop_duplicates(['TSS'])
	sub_gtf_df['TSS'] = sub_gtf_df['TSS'].astype(int)

	real_ep_df = ep_df[['chrom-Enh','chromStart','chromEnd','TSS']]
	real_ep_df = real_ep_df.sort_values(by=['chromEnd'],ascending=True)#.drop_duplicates(['chromStart','chromEnd'])
	real_ep_df[['chromStart','chromEnd']] = real_ep_df[['chromStart','chromEnd']].astype(int)

	arrived_index = 0;
	sub_gtf_df_tss_loc_vec = sub_gtf_df['TSS'].values
	fake_ep_df = []
	for item in real_ep_df.values:
		tmp_fake_ep = [item[0],item[1],item[2]] 
		for gtf_tss in sub_gtf_df_tss_loc_vec[arrived_index:]:
			right_enhancer_distance = abs(item[2] - gtf_tss)
			left_enhancer_distance = abs(item[1]-gtf_tss)
			if (right_enhancer_distance>4000 and right_enhancer_distance<20000) or\
			   (left_enhancer_distance>4000 and left_enhancer_distance<20000):
				tmp_fake_ep.append(gtf_tss)
				fake_ep_df.append(tmp_fake_ep)
				break
			arrived_index+=1
	fake_ep_df = pd.DataFrame(fake_ep_df,columns=['chrom-Enh','chromStart','chromEnd','TSS'])
	print(fake_ep_df.shape)
	print(real_ep_df.shape)

	fake_ep_df['label'] = 0
	real_ep_df['label'] = 1
	real_feak_ep_df = real_ep_df.append(fake_ep_df)
	print(real_feak_ep_df.shape)
	print(real_feak_ep_df.isna().any())
	real_feak_ep_df.to_csv(save_path,index=False)    
	print('cell {} fake ep and real ep has been matched'.format(celline))

if __name__ == '__main__':
	generate_fake_ep('K562')