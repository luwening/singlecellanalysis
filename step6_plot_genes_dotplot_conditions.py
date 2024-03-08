#usage: python3 ~/pipeline/single_cell/10x/scanpy_pipeline/step7_plot_genes_volin.py step2_scanpy_temp_file_merged_14_samples_bbknn.h5ad fibroblast_genes
import numpy as np
import pandas as pd 
import scanpy as sc 
import matplotlib.pylab as plt 
#import doubletdetection
import os

sc.settings.verbosity=3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80,facecolor='white')
sc.settings.n_jobs=10

import sys

name=sys.argv[1] # input name
gene_name=sys.argv[2]

adata=sc.read(name)

#marker_genes=['DCN','GSN','APOD','CXCL12','CLDN4','TACSTD2','TM4SF1','MGP','NKG7','B2M','HLA-B','HLA-A','VIM', 'CRYAB','S100B','LGALS1']
marker_genes=list(set(open(gene_name).read().split()))
marker_genes=open(gene_name).read().split()

#vp=sc.pl.dotplot(adata,marker_genes,groupby='sample',return_fig=True,dendrogram=True)
vp=sc.pl.dotplot(adata,marker_genes,groupby='condition',return_fig=True,figsize=(8,2.5))

#vp.add_totals().style(ylim=(0,8))
vp.savefig('step7_%s_%s_marker_genes_dotplot_condition.png' % (gene_name,os.path.basename(name)),dpi=300)