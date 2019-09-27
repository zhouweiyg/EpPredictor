import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import accuracy_score,f1_score
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Imputer
from sklearn.externals import joblib
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

import numpy as np

import lightgbm as lgb

DROP_COL = ['chromStart','chromEnd','id','chrom-Enh','TSS','label']


def lable_hic(features, hi_c_train):
	print(features.shape, hi_c_train.shape)
	scaler = MinMaxScaler()
	distance = hi_c_train['hic_distance'].values
	distance = np.reshape(distance, (-1, 1))
	label = scaler.fit_transform(distance)
	label = label[:, 0]
	hi_c_train['hic_distance'] = label
	hi_c_train = hi_c_train.rename(columns={'hi_c_train':'label'})
	features = features.drop('label', axis=1)
	features = features.merge(hi_c_train, how='left', on=['chrom-Enh','chromStart','chromEnd','TSS'])
	
	return features
	

def load_data(features):
	"""
	load complete feature and lables, then split the data to test and train

	@features: complete features in pandas DataFrame format
	@label: 5D vector, the first column is its id

	return: train and test data set while the label which I used it mean value
	"""

	y = features['label']
	X = features.drop(columns=DROP_COL, axis=1)
	X_train, X_test, y_train, y_test = train_test_split(
		X, y, test_size=0.33, random_state=42)
    
	return X_train, X_test, y_train, y_test



if __name__ == '__main__':
	features = pd.read_csv('../midfile/train2.csv')
	hi_c_train = pd.read_csv('../midfile/hic-train.csv')
	hi_c_train = hi_c_train.fillna(0)

	features = features.fillna(method='pad')

	lable_hic(features, hi_c_train)