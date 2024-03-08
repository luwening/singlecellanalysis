# merge the concat data using bbknn method
# usage: python3 ~/pipeline/single_cell/10x/scanpy_pipeline/step1_merge_bbknn.py

import numpy as np
import pandas as pd 
import scanpy as sc 
import matplotlib.pylab as plt 
import os
import bbknn

#import doubletdetection

sc.settings.verbosity=3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80,facecolor='white')
sc.settings.n_jobs=10

import sys

name=sys.argv[1]

#adata_concat=sc.read('scanpy_temp_file_merged_14_samples_simple_concat.h5ad')
adata_concat=sc.read(name)
adata_concat.layers['counts']=adata_concat.X.copy()

sc.pp.normalize_total(adata_concat,target_sum=1e4)
sc.pp.log1p(adata_concat)

sc.tl.pca(adata_concat)
sc.pp.neighbors(adata_concat)
sc.tl.umap(adata_concat)

sc.set_figure_params(figsize=(8, 5))

fig=sc.pl.umap(adata_concat, color=[ 'celltype'],ncols=1,add_outline=True,return_fig=True,palette=sc.pl.palettes.vega_20_scanpy,legend_fontsize=12,legend_fontoutline=2,frameon=True,title='clustring of cells')
#sc.set_figure_params(dpi=80,figsize=(8,5))
#more right margin
#fig.subplots_adjust(left=0, bottom=0, right=0.5, top=1, wspace=0, hspace=0)
plt.tight_layout(pad=0.05)

plt.savefig('%s_umap.png' % os.path.basename(name).split('.h5ad')[0],dpi=300)

adata_concat.layers['log1norm']=adata_concat.X
adata_concat.X=adata_concat.X.toarray()

print('max,min,median value in adata_concat:')
print(np.max(adata_concat.X),np.min(adata_concat.X),np.median(adata_concat.X))

bbknn.bbknn(adata_concat, batch_key='sample')
print('max,min,median value in adata_concat:')
print(np.max(adata_concat.X),np.min(adata_concat.X),np.median(adata_concat.X))
sc.tl.umap(adata_concat)
sc.set_figure_params(figsize=(6,5))
sc.pl.umap(adata_concat, color=['sample'], palette=sc.pl.palettes.vega_20_scanpy,ncols=1)
#sc.set_figure_params(dpi=80,figsize=(6,5))
plt.tight_layout(pad=0.05)
plt.savefig('%s_bbknn_umap_sample.png' % os.path.basename(name).split('.h5ad')[0])

sc.set_figure_params(figsize=(6,5))
sc.pl.umap(adata_concat, color=[ 'celltype'], palette=sc.pl.palettes.vega_20_scanpy,ncols=1)
#sc.set_figure_params(dpi=80,figsize=(6,5))
plt.tight_layout(pad=0.05)
plt.savefig('%s_bbknn_umap_celltype.png' % os.path.basename(name).split('.h5ad')[0])

sc.tl.leiden(adata_concat,  key_added="leiden_res0_25", resolution=0.25)
sc.tl.leiden(adata_concat,  key_added="leiden_res0_4", resolution=0.4)
sc.tl.leiden(adata_concat,  key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata_concat,  key_added="leiden_res1", resolution=1)
sc.pl.umap(adata_concat, color = ["leiden_res0_25","leiden_res0_4","leiden_res0_5","leiden_res1"],legend_loc='on data', palette=sc.pl.palettes.vega_20_scanpy,ncols=1)
plt.savefig('%s_0.4_leiden.png' % os.path.basename(name).split('.h5ad')[0])

#convert matrix to np.array
# max,min,median value in adata_concat
print('max,min,median value in adata_concat:')
print(np.max(adata_concat.X),np.min(adata_concat.X),np.median(adata_concat.X))
bbknn.ridge_regression(adata_concat, batch_key=['sample'], confounder_key=['leiden_res0_4']) # need or not?
adata_concat.layers['ridge_data']=adata_concat.X
print('max,min,median value in adata_concat:')
print(np.max(adata_concat.X),np.min(adata_concat.X),np.median(adata_concat.X))
sc.pp.pca(adata_concat)

bbknn.bbknn(adata_concat, batch_key='sample')
print('max,min,median value in adata_concat:')
print(np.max(adata_concat.X),np.min(adata_concat.X),np.median(adata_concat.X))
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['sample','celltype'],ncols=1, palette=sc.pl.palettes.vega_20_scanpy)
#sc.set_figure_params(dpi=80,figsize=(8,5))
plt.savefig('%s_bbknn_umap_ridge.png' % os.path.basename(name).split('.h5ad')[0])


adata_concat.X=adata_concat.layers['counts']
result_file='%s_bbknn.h5ad' % os.path.basename(name).split('.h5ad')[0]
adata_concat.write(result_file)

name='celltype_name_bbknn_reduce_inhouse_preprocessed_nonorm_11samples.h5ad'

result_file='%s_bbknn.h5ad' % os.path.basename(name).split('.h5ad')[0]

adata_concat=sc.read(result_file)

adata_concat.X=adata_concat.layers['log1norm']
#highly variable genes
if 'base' not in adata_concat.uns['log1p']:
    adata_concat.uns['log1p']['base'] = None
sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_concat)
plt.savefig('%s_highly_variable_genes.png' % os.path.basename(name).split('.h5ad')[0],dpi=100)

sc.tl.rank_genes_groups(adata_concat,'leiden_res0_4',method='wilcoxon')

sc.pl.rank_genes_groups(adata_concat,n_genes=25,sharey=False)
plt.savefig('%s_rank_genes.png' % os.path.basename(name).split('.h5ad')[0],dpi=100)


# add condition
# if sample in 0,2,4,6,8,10,12, tumor; else normal
adata_concat.obs['condition'] = 'normal'
adata_concat.obs.loc[adata_concat.obs['sample'].isin(['0','2','4','6','8','10','12']), 'condition'] = 'tumor'


marker_genes=['DCN','CFD','CD3D']


marker_genes_dict={
    'cancercells':['EPCAM','KRT19','KRT14'],
    'cancerstemcells':['KRT19','TOP2A'],
    'CXCL14cancer':['CXCL14'],
    'Bcells':['CD79A','CD79B'],
    'Macrophages':['LYZ','IL1B','MSR1'],
    'myeloid':['LYZ'],
    'matrueDC':['LAMP3','CCR7'],
    'plasma':['JCHAIN','MZB1'],
    'endothelial':['MCAM','PECAM1'],
    'pericytes':['ACTA2','TAGLN','MCAM'],
    'myofibroblasts':['DCN','TAGLN'],
    'cyclingcells':['TOP2A','MKI67'],
    'epithelialcells':['VWF','MCAM','ERBB2','ESR1'],
    'naiveT':['CD3G','CD3D','IL7R','NKG7'],
    'CD8effector':['CD3D','IL7R','NKG7','GNLY','CD8A']
    }
marker_genes_dict={
    'cancercells':['EPCAM','KRT19','KRT14','TOP2A'],
    'Bcells':['CD79A','CD79B'],
    'Macrophages':['LYZ','IL1B','MSR1'],
    'myeloid':['LYZ'],
    'matrueDC':['LAMP3','CCR7'],
    'plasma':['JCHAIN','MZB1'],
    'endothelial':['MCAM','PECAM1'],
    'pericytes':['ACTA2','TAGLN','MCAM'],
    'myofibroblasts':['DCN','TAGLN'],
    'epithelialcells':['VWF','MCAM','ERBB2','ESR1'],
    'Tcells':['CD3G','CD3D','IL7R','NKG7','GNLY','CD8A'],
    }
sc.pl.umap(adata_concat, color=["condition", "celltype"], frameon=False)
plt.savefig('%s_condition.png' % os.path.basename(name).split('.h5ad')[0],dpi=100)

output_name_prefix=os.path.basename(name).split('.h5ad')[0]

vp=sc.pl.stacked_violin(adata_concat,marker_genes,groupby='leiden_res0_4',rotation=90,return_fig=True)
vp.add_totals().style(ylim=(0,5))
#vp.savefig('%s_marker_genes_volin.png' % sample,dpi=300)
vp.savefig('%s_rank_genes_wilcoxon.png' % output_name_prefix,dpi=100)

sc.pl.rank_genes_groups_stacked_violin(adata_concat,n_genes=6,cmap='viridis_r')
plt.tight_layout()

plt.savefig('%s_stacked_violin.png' % output_name_prefix,dpi=200)


rank_genes_groups=adata_concat.uns['rank_genes_groups']
groups=rank_genes_groups['names'].dtype.names
marker=pd.DataFrame({group+'_'+key[:1]: rank_genes_groups[key][group] for group in groups for key in ['names','pvals']})
marker.to_csv('%s_marker_leiden.csv' % output_name_prefix)


marker_genes_dict={
    'cancercells':['EPCAM','KRT19','KRT14','TOP2A','VWF','MCAM','ERBB2','ESR1'],
    'Bcells':['CD79A','CD79B'],
    'Macrophages':['LYZ','IL1B','MSR1'],
    'plasma':['JCHAIN','MZB1'],
    'endothelial':['MCAM','PECAM1'],
    'pericytes':['ACTA2','TAGLN','MCAM'],
    'myofibroblasts':['DCN','TAGLN','COL1A1','COL1A2','COL3A1','CFD'],
    'Tcells':['CD3G','CD3D','IL7R','NKG7','GNLY','CD8A'],
    }

sc.pl.dotplot(adata_concat,marker_genes_dict,'leiden_res0_4',dendrogram=True,figsize=(16,8))
plt.tight_layout(pad=10,h_pad=8)
#fig.subplots_adjust(left=0.2, right=0.98, top=0.9, bottom=0.1)
plt.savefig('%s_marker_gene_dict_dotplot.png' % output_name_prefix)




cluster2annotation={
    #total 41: 0-40
    '1':'Fibroblasts',

    '5':'Perivascular cells', # not sure


    '8':'B cells',



    '10':'Macrophages',


    '2':'Epithelial cells', #important?
    '3':'Epithelial cells', #important?
    '9':'Epithelial cells', #important?
    '6':'Epithelial cells', # not sure


    '0':'T cells',
    '7':'T cells',


    '11':'Plasma cells',





    '4':'Endothelial cells',
 

    #'22':'myofibroblasts',
    #'23':'',
}

names_new={'Tcells':'T cells','myofibroblasts':'Fibroblasts','cancercells':'Epithelial cells','endothelial':'Endothelial cells','pericytes':'Pericytes','plasma':'Plasma cells','Macrophages':'Macrophages','Bcells':'B cells'}


for temp_item in range(12):
    print(temp_item,cluster2annotation[str(temp_item)])

adata_concat.obs['celltype_new']=adata_concat.obs['leiden_res0_4'].map(cluster2annotation).astype('category')
sc.pl.umap(adata_concat,color='celltype_new',legend_loc='on data',frameon=False,legend_fontsize=12,legend_fontoutline=2,title='cell type',save='%s_leiden_celltype.png' % name)

adata_concat.write('%s_population.h5ad' % name)

#paga
'''
adata_new=adata_concat[(adata_concat.obs['celltype']=='CXCL14cancer')|(adata_concat.obs['celltype']=='cancercells')|(adata_concat.obs['celltype']=='cancerstemcells'),:]
sc.tl.paga(adata_new,groups='celltype')
sc.pl.paga(adata_new,color='celltype',save='cancer_cell_paga.png')
sc.tl.draw_graph(adata_new,init_pos='paga')
adata_new.uns['iroot']=np.flatnonzero(adata_new.obs['celltype']=='cancerstemcells')[0]
sc.tl.dpt(adata_new)
sc.pl.draw_graph(adata_new,color=['celltype','dpt_pseudotime','group'],legend_loc='on data',save='cancercells_faga_new.png')

sc.tl.paga(adata,groups='celltype')

sc.pl.paga(adata,color=['celltype'],node_size_scale=0.5,save='cancercells_faga_new2.png')
sc.tl.draw_graph(adata,init_pos='paga')
adata.uns['iroot']=np.flatnonzero(adata.obs['celltype']=='NaiveT')[0]
sc.tl.dpt(adata)
sc.pl.draw_graph(adata,color=['celltype','dpt_pseudotime','group'],legend_loc='on data',save='all_cells_faga_new.png')
'''


