import numpy as numpy
import pandas as pd 
import scanpy as sc 
import matplotlib.pylab as plt 
import os
import sys
import numpy as np
#import doubletdetection
#import pandas

#pandas.set_option('mode.use_inf_as_na', True)


sc.settings.verbosity=3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80,facecolor='white')
sc.settings.n_jobs=10


#prefix='simple_concat'
#prefix='combat'
#prefix='bbknn'

name=sys.argv[1]
basename=os.path.basename(name).split('.')[0]


#result_file='scanpy_temp_file_merged_14_samples_%s.h5ad' % prefix
result_file=name
#~/pipeline/single_cell/10x/scanpy_pipeline/step2_populations.py
adata=sc.read(name)

#sample='merged_14_samples_%s' % prefix
sample=basename
#adata=adata_merge
#adata.X.data=np.nan_to_num(adata.X.data)
if 'log1p' in adata.uns:
    adata.uns['log1p'] ["base"] = None
#adata.var['mt']=adata.var_names.str.startswith('MT-')
#sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
sc.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=3,min_disp=0.5)
sc.pl.highly_variable_genes(adata)
plt.savefig('%s_highly_variable_genes.png' % sample)

adata.raw=adata
adata=adata[:,adata.var.highly_variable]


#sc.pp.regress_out(adata,['total_counts','pct_counts_mt']) # need or not?
#sc.pp.scale(adata,max_value=10)

sc.tl.pca(adata,svd_solver='arpack')
sc.pl.pca(adata,color='DCN')
plt.savefig('%s_PCA.png' % sample)

sc.pl.pca_variance_ratio(adata,log=True)
plt.savefig('%s_variance_ratio.png' % sample)


marker_genes=['DCN','CFD']


sc.pp.neighbors(adata,n_neighbors=20,n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata,color=marker_genes)
plt.savefig('%s_umap.png' % sample)

sc.tl.leiden(adata,key_added='clusters',resolution=0.5)
sc.pl.umap(adata,color='clusters',add_outline=True,legend_loc='on data',use_raw=False, legend_fontsize=12,legend_fontoutline=2,frameon=False,title='clustering of cells',palette='Set1')
plt.savefig('%s_umap_leiden.png' % sample ,dpi=200)


#sc.tl.rank_genes_groups(adata,'clusters',methods='t-test')  #different test for rank genes
sc.tl.rank_genes_groups(adata,'clusters',methods='wilcoxon')

sc.pl.rank_genes_groups(adata,n_genes=25,sharey=False)
plt.savefig('%s_rank_genes.png' % sample,dpi=100)
vp=sc.pl.stacked_violin(adata,marker_genes,groupby='clusters',rotation=90,return_fig=True)
vp.add_totals().style(ylim=(0,5))
#vp.savefig('%s_marker_genes_volin.png' % sample,dpi=300)
vp.savefig('%s_rank_genes_wilcoxon.png' % sample,dpi=100)

sc.pl.rank_genes_groups_stacked_violin(adata,n_genes=6,cmap='viridis_r')
plt.savefig('%s_stacked_violin.png' % sample,dpi=200)

#result_file='step2_scanpy_temp_file_%s.h5ad' % sample
adata.write('%s_%s.h5ad' %(basename,'population'))