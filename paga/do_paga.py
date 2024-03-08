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

adata=sc.read(name)

#adata_new=adata[(adata.obs['celltype']=='CXCL14cancer')|(adata.obs['celltype']=='cancercells')|(adata.obs['celltype']=='cancerstemcells'),:]
adata_new=adata

sc.tl.paga(adata_new,groups='celltype')
sc.pl.paga(adata_new,color='celltype',save='cancer_cell_paga.png')
sc.tl.draw_graph(adata_new,init_pos='paga')
adata_new.uns['iroot']=np.flatnonzero(adata_new.obs['celltype']=='cancercells')[1000]
sc.tl.dpt(adata_new)
sc.pl.draw_graph(adata_new,color=['celltype','dpt_pseudotime','group'],legend_loc='on data',save='cancercells_faga_new.png')

sc.tl.paga(adata,groups='celltype')

sc.pl.paga(adata,color=['celltype'],node_size_scale=0.5,save='cancercells_faga_new2.png')
sc.tl.draw_graph(adata,init_pos='paga')
adata.uns['iroot']=np.flatnonzero(adata.obs['celltype']=='NaiveT')[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata,color=['celltype','dpt_pseudotime','group'],legend_loc='on data',save='all_cells_faga_new.png')
