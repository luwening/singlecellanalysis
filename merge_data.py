import sys
import scanpy as sc
import numpy as np

#working dir: /disk3/users/lvwen/download/inhouse/merge_3datasets
names='''
/disk3/users/lvwen/download/inhouse/bbknn/scanpy_temp_file_merged_14_samples_simple_concat_nopreprocess.h5ad
/disk3/users/lvwen/download/GSE164898/scanpy_temp_file_merged_14_samples_simple_concat.h5ad
/disk3/users/lvwen/download/GSE180286/data_merge.h5ad'''.split()

data={}
for name in names:
    data[name]=sc.read(name)

samples=list(data.keys())

adata_merge=data[samples[0]].concatenate(data[samples[1]],data[samples[2]],join='outer')
adata_merge.X.data=np.nan_to_num(adata_merge.X.data)

#sc.pp.filter_genes(adata_merge,min_cells=5) # filter gene here?
#result_file='scanpy_temp_file_%s_simple_concat.h5ad' % 'merged_14_samples'
result_file='name_bbknn.h5ad'

#data_merge.var=adata_merge.var.drop(columns=['mt-1'])
adata_merge.write(result_file)