# -*- coding: utf-8 -*-
"""
scaler for matrix with header and index
created on 2023-04-01 @Guo
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import QuantileTransformer
import sys,re

data=pd.read_csv(sys.argv[1],sep='\t',header=0,index_col=0)
header=data.columns
index=data.index

data=np.array(data)

scaler=MinMaxScaler()
scaler.fit(data)
X=pd.DataFrame(scaler.transform(data))
X.index=index
X.columns=header

qt=QuantileTransformer(n_quantiles=1000, random_state=0)
Y=pd.DataFrame(qt.fit_transform(data))
Y.index=index
Y.columns=header

X.to_csv(sys.argv[2], sep='\t')
Y.to_csv(sys.argv[3], sep='\t')
