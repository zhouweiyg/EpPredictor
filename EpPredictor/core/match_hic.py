import pandas as pd
import numpy as np

EP_COLUMNS_NAME = ['chrom-Enh','chromStart','chromEnd','Gene ID','chrom-Gene','TSS','Transcript',
					'signalValue','EP Score']

EP_COLUMNS_NAME_NEW = ['chrom-Enh','chromStart','chromEnd','Gene ID','chrom-Gene','TSS','Transcript',
					'signalValue','EP Score','label']

TAD_COLUMNS_NAME = ['chrom','tad_start','tad_end']


def load_ep(ep_fpath, seps):
	ep_df = pd.read_csv(ep_fpath, sep=seps, names=EP_COLUMNS_NAME)
	return ep_df

def load_tad(tad_fpath):
	tad_df = pd.read_csv(tad_fpath, sep='\t', names=TAD_COLUMNS_NAME)
	return tad_df


def process(ep_df):
	"""
	chr start end chr start end color commnet
	"""
	put_ep = pd.DataFrame()

	put_ep['chr1'] = ep_df['chrom-Enh']
	put_ep['x1'] = (ep_df['chromEnd']+ep_df['chromStart'])//2
	put_ep['x2'] = ep_df['TSS']
	put_ep['chr2'] = ep_df['chrom-Enh']
	put_ep['y1'] = put_ep['x1']
	put_ep['y2'] = put_ep['x2']
	put_ep['color'] = '0,255,0'
	put_ep['comment'] = 'my gren region'
	
	return put_ep


def process_tad(tad_df):
	"""
	"""
	put_tad = pd.DataFrame()

	put_tad['chr1'] = tad_df['chrom']
	put_tad['x1'] = tad_df['tad_start']
	put_tad['x2'] = tad_df['tad_end']
	put_tad['chr2'] = tad_df['chrom']
	put_tad['y1'] = put_tad['x1']
	put_tad['y2'] = put_tad['x1']
	put_tad['color'] = '0,0,255'
	put_tad['comment'] = 'my gblue region'

	return put_tad


def match_label(ep_df, hic_map_df, ep_hic_save):
	hic_map_df['chr1'] = 'chr'+hic_map_df['chr1']

	std_ep_lable = []

	for ep in ep_df.values:
		ep_chrom = ep[0]
		ep_x1 = int((ep[1]+ep[2])//2)
		ep_x2 = int(ep[5])
		if ep_x1>ep_x2:
			ep_x1, ep_x2 =ep_x2, ep_x1
		lable = 0

		current_hic_map_df = hic_map_df[hic_map_df['chr1']==ep_chrom]
		for hic_item in current_hic_map_df.values:
			hic_x1 = int(hic_item[1])
			hic_x2 = int(hic_item[2])

			if(ep_x1>=hic_x1 and ep_x2<=hic_x2):
				lable = 1
				break;
		ep = np.r_[ep, lable]
		std_ep_lable.append(ep)

	std_ep_lable = pd.DataFrame(std_ep_lable, columns=EP_COLUMNS_NAME_NEW)
	std_ep_lable.to_csv(ep_hic_save,index=False)

'''

if __name__ == '__main__':
	ep_df = load_ep('./data/K562_EP.txt')
	
	#put_ep = process(ep_df)
	#put_ep.to_csv('./midfile/put_ep.csv', index=False, sep='\t')

	#tad_df = load_tad('./data/K562_Lieberman-raw_TADs.txt')
	#put_tad = process_tad(tad_df)
	#put_tad.to_csv('./midfile/put_tad.csv', index=False, sep='\t')

	hic_map_df = pd.read_csv('./data/hic-map.txt', sep='\t')
	match_label(ep_df, hic_map_df)
'''
