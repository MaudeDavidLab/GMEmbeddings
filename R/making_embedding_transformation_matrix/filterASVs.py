#!/usr/bin/env python

import pandas as pd
import time
import numpy as np

chunks = pd.read_csv('/nfs3/PHARM/David_Lab/austin_e/filtered_final_seqtab/final_seqtab2_1.csv', chunksize=100)

start = time.time()
df = pd.DataFrame([])
for chunk in chunks:
    temp_df = chunk
    print("Initial size of chunk: ", temp_df.shape)
    temp_df.replace(0, np.nan, inplace=True) 
    num_count = temp_df.count(axis=1) 
    filtered_index = num_count[num_count > 10].index 
    temp_df.loc[temp_df.index.isin(filtered_index)]
    output_df = temp_df.loc[temp_df.index.isin(filtered_index)]
    output_df.fillna(0, inplace=True) 
    print("Final size of chunk: ", output_df.shape)
    df = pd.concat([df, output_df])
    print("Size of concatenated df: ", df.shape)
end = time.time()
print("total time to complete: ",(end-start), "sec")

df.to_csv(r'/nfs3/PHARM/David_Lab/austin_e/filtered_final_seqtab/filtered_final_seqtab2.csv')
