import numpy as np
import pandas as pd 
import scanpy as sc 
import matplotlib.pylab as plt 
import glob
import os 
import sys

sc.settings.verbosity=3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80,facecolor='white')
sc.settings.n_jobs=10


name=sys.argv[1]

adata=sc.read_h5ad(name)
adata.X.data=np.nan_to_num(adata.X.data)
#print(adata)
#print(adata.X)

#need to think about how to deal with the missing  value
#adata.X=np.nan_to_num(adata.X)
#adata=adata.transpose()

#adata=sc.read_10x_h5(name)
adata.var_names_make_unique()

print(adata)
sc.pl.highest_expr_genes(adata,n_top=50)
sample=name.split('.')[0]
plt.tight_layout()
plt.savefig('%s_high_expressed_genes.png' % sample)

#do not filter before merge data
sc.pp.filter_cells(adata,min_genes=200)
sc.pp.filter_genes(adata,min_cells=3)

print(adata)
#print(adata.X)

adata.var['mt']=adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
sc.pl.violin(adata,['n_genes_by_counts','total_counts','pct_counts_mt'],jitter=0.4,multi_panel=True)
plt.savefig('%s_mt_percentage_before_filter.png' % sample)

sc.pl.scatter(adata,x='total_counts',y='pct_counts_mt')
plt.savefig('%s_counts_vs_mt.png' % name)

sc.pl.scatter(adata,x='total_counts',y='n_genes_by_counts')
plt.savefig('%s_counts_vs_n_genes.png' % name)

adata=adata[adata.obs.pct_counts_mt<20,:]
adata=adata[adata.obs.n_genes_by_counts<5000,:]

sc.pl.violin(adata,['n_genes_by_counts','total_counts','pct_counts_mt'],jitter=0.4,multi_panel=True)
plt.savefig('%s_mt_percentage.png' % sample)



groups=[]
for cell in adata.obs.index:
    if '_' in cell:
        groups.append(cell.split('_',1)[1])
    elif '-' in cell:
        groups.append(cell.split('-',1)[1])

for group in set(groups):
    print(group,groups.count(group))

adata.obs['group']=groups

adata.write('%s_preprocessed_nonorm.h5ad' % sample)

sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)


adata.write('%s_preprocessed.h5ad' % sample)