import pandas as pd
import sys
import os

store =pd.HDFStore('../results/final_table.HDF', complib='blosc')


for file in os.listdir('../tmp/results/'):
    df=pd.read_hdf('../tmp/results/'+file)
    store.append(key='data',value=df,format='t',data_columns=['Homo_sapiens'],
                 min_itemsize={ 'chr.Homo_sapiens' : 30 ,
                                'chr.Gorilla_gorilla_gorilla': 30,
                                'chr.Pan_troglodytes': 30,
                                'chr.Pongo_abelii': 30})

store.close()
