#usage: python3 /home/lu/pipeline/single_cell/10x/scanpy_pipeline/step7_plot_genes_umap.py name_bbknn_reduce_inhouse_preprocessed_population.h5ad breast_cancer_celltype_markers

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
import matplotlib as mpl
from copy import copy

reds = copy(mpl.cm.Reds)

name=sys.argv[1] # input name
gene_name=sys.argv[2]

adata=sc.read(name)


if 'log1p' in adata.layers.keys():
    adata.X=adata.layers['log1p']
#marker_genes=['DCN','GSN','APOD','CXCL12','CLDN4','TACSTD2','TM4SF1','MGP','NKG7','B2M','HLA-B','HLA-A','VIM', 'CRYAB','S100B','LGALS1']
marker_genes=open(gene_name).read().split()

new_genes=[]
for gene in marker_genes:
    if gene in adata.var.index:
        new_genes.append(gene)

print(set(marker_genes)-set(new_genes))
#vp=sc.pl.stacked_violin(adata,marker_genes,groupby='clusters',rotation=90,return_fig=True,dendrogram=True)
#vp=sc.pl.umap(adata,color=new_genes,return_fig=True)

# cmap gray to red
cmap=mpl.colors.LinearSegmentedColormap.from_list('mycmap', ['lightgreen', 'red'])
vp=sc.pl.umap(adata,color=marker_genes,return_fig=True,ncols=3,cmap=cmap,vmin=0, add_outline=True,frameon=False)

#plt.savefig('%s_umap.png' % sample)
#vp.add_totals().style(ylim=(0,8))
vp.savefig('step7_%s_%s_marker_genes_umap.png' % (os.path.basename(gene_name),os.path.basename(name)),dpi=300)