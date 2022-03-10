import pandas as pd
import numpy as np
from sklearn import preprocessing

data = pd.read_csv("seqtab_match_mapping.csv", index_col = 0, header = 0)
#data = data.iloc[0:10, 0:10]
data_center = preprocessing.scale(data)
u,s,vh = np.linalg.svd(data_center, full_matrices = False)

v = np.transpose(vh)
v = pd.DataFrame(v)
v.index = data.columns.values
v.to_csv("pca_transformation_matrix.csv")

