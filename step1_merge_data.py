import numpy as np
import pandas as pd 
import scanpy as sc 
import matplotlib.pylab as plt 

import doubletdetection

sc.settings.verbosity=3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80,facecolor='white')
sc.settings.n_jobs=10


#read data
def process_sample_normalized(sample):
    adata=sc.read_10x_mtx('/data/single_cell_data/in_lab_7_pair/%s/' % sample,var_names='gene_symbols',cache=True)
    adata.var_names_make_unique()

    print(adata)
    #sc.pl.highest_expr_genes(adata,n_top=30)
    #plt.savefig('%s_high_expressed_genes.png' % sample)

    #do not filter before merge data
    sc.pp.filter_cells(adata,min_genes=200)
    #sc.pp.filter_genes(adata,min_cells=3)
    adata.var['mt']=adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)

    adata=adata[adata.obs.pct_counts_mt<20,:]
    adata=adata[adata.obs.n_genes_by_counts<8000,:]

    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)

    #sc.pl.violin(adata,['n_genes_by_counts','total_counts','pct_counts_mt'],jitter=0.4,multi_panel=True)
    #plt.savefig('%s_mt_percentage.png' % sample)
    return adata

samples=['X1Ca','X1Pa','X2Ca','X2Pa','X3Ca','X3Pa','X4Ca','X4Pa','X5Ca','X5Pa','X6Ca','X6Pa','X7Ca','X7Pa']
target_samples=['X1Ca','X1Pa','X2Ca','X2Pa','X4Ca','X4Pa','X5Ca','X6Ca','X6Pa','X7Ca','X7Pa']


all_data={}
for sample in samples:
    all_data[sample]=process_sample_normalized(sample)
    #check_doublet(sample)

adata_merge=all_data[samples[0]].concatenate(all_data[samples[1]],all_data[samples[2]],all_data[samples[3]],all_data[samples[4]],all_data[samples[5]],all_data[samples[6]],all_data[samples[7]],all_data[samples[8]],all_data[samples[9]],all_data[samples[10]],all_data[samples[11]],all_data[samples[12]],all_data[samples[13]],batch_key='sample')
sc.pp.filter_genes(adata_merge,min_cells=5) # filter gene here?
result_file='scanpy_temp_file_%s_simple_concat_nopreprocess.h5ad' % 'merged_14_samples'
adata_merge.write(result_file)

sc.pp.combat(adata_merge,key='sample')
result_file='scanpy_temp_file_%s_nopreprocess.h5ad' % 'merged_14_samples'
adata_merge.write(result_file)