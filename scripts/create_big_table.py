import pandas as pd
import sys
import os
import pickle

target_seqname = sys.argv[1]

slice_lst = pickle.load(open('../tmp/slice_dct.txt', 'rb'))

store = pd.HDFStore('../results/final_table.HDF', complib='blosc')

for run in range(len(slice_lst)):
    df = pd.read_hdf('../tmp/results/run_{}.HDF'.format(run))
    store.append(key=target_seqname,value=df,format='t',data_columns=['Homo_sapiens'],
                 min_itemsize={ 'chr.Homo_sapiens' : 30 ,
                                'chr.Gorilla_gorilla_gorilla': 30,
                                'chr.Pan_troglodytes': 30,
                                'chr.Pongo_abelii': 30})

store.close()
