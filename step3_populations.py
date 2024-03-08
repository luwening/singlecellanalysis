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
sample=name

adata=sc.read_h5ad(name)
adata.uns['log1p']["base"] = None
sc.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=3,min_disp=0.5)
sc.pl.highly_variable_genes(adata)
plt.savefig('%s_highly_variable_genes.png' % sample)



adata.raw=adata
adata=adata[:,adata.var.highly_variable]
#sc.pp.regress_out(adata,['total_counts','pct_counts_mt']) # need or not?
#sc.pp.scale(adata,max_value=10)

print(adata.shape)

#still need?
#sc.pp.scale(adata,max_value=10)



sc.tl.pca(adata,svd_solver='arpack')
sc.pl.pca(adata,color='CD3D')
plt.savefig('%s_PCA.png' % sample)

sc.pl.pca_variance_ratio(adata,log=True)
plt.savefig('%s_variance_ratio.png' % sample)


marker_genes=['DCN','CFD','CD3D']


sc.pp.neighbors(adata,n_neighbors=20,n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata,color=marker_genes)
plt.savefig('%s_marker_gene_umap.png' % sample)

sc.tl.umap(adata)
#sc.pl.umap(adata,color='group')

with plt.rc_context({'figure.figsize':(6,6)}):
    sc.pl.umap(adata,color='group',add_outline=True,legend_fontsize=12,legend_fontoutline=2,frameon=False,title='clustring of cells',palette='Paired')
    plt.tight_layout()
    plt.savefig('%s_group_umap.png' % sample)


sc.tl.leiden(adata)
sc.tl.louvain(adata)

sc.pl.umap(adata,color=['leiden','louvain'],add_outline=True,legend_fontsize=5,legend_fontoutline=2,frameon=False,title=['leiden','louvain'],palette='Set1',legend_loc='on data')
plt.tight_layout()
plt.savefig('%s_clusters_umap.png' % sample)

sc.tl.rank_genes_groups(adata,'leiden',method='wilcoxon')

sc.pl.rank_genes_groups(adata,n_genes=25,sharey=False)
plt.savefig('%s_rank_genes.png' % sample,dpi=100)

vp=sc.pl.stacked_violin(adata,marker_genes,groupby='leiden',rotation=90,return_fig=True)
vp.add_totals().style(ylim=(0,5))
#vp.savefig('%s_marker_genes_volin.png' % sample,dpi=300)
vp.savefig('%s_rank_genes_wilcoxon.png' % sample,dpi=100)

sc.pl.rank_genes_groups_stacked_violin(adata,n_genes=6,cmap='viridis_r')
plt.tight_layout()

plt.savefig('%s_stacked_violin.png' % sample,dpi=200)


rank_genes_groups=adata.uns['rank_genes_groups']
groups=rank_genes_groups['names'].dtype.names
marker=pd.DataFrame({group+'_'+key[:1]: rank_genes_groups[key][group] for group in groups for key in ['names','pvals']})
marker.to_csv('%s_marker_leiden.csv' % name)

marker_genes_dict={
    'cancercells':['EPCAM','KRT19','KRT14'],
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

marker_genes_dict2={
    'cancercells':['EPCAM','KRT19','KRT14'],
    'cancerstemcells':['KRT19','TOP2A'],
    'CXCL14cancer':['CXCL14'],
    'Bcells':['CD79A','CD79B'],
    'Macrophages':['LYZ','IL1B','MSR1'],
    'myeloid':['LYZ'],
    'matrueDC':['LAMP3','CCR7'],
    'plasma':['JCHAIN','IGHG3','MZB1'],
    'endothelial':['MCAM','PECAM1'],
    'pericytes':['ACTA2','TAGLN','MCAM'],
    'myofibroblasts':['DCN','LUM','TAGLN'],
    'cyclingcells':['TOP2A','MKI67'],
    'epithelialcells':['VWF','MCAM','ERBB2','ESR1'],
    'naiveT':['CD3G','CD3D','IL7R','NKG7'],
    'CD8effector':['CD3D','IL7R','NKG7','GNLY','CD8A']
    }

sc.pl.stacked_violin(adata,marker_genes_dict,groupby='leiden',swap_axes=False,dendrogram=True,cmap='Paired_r')
plt.tight_layout()

plt.savefig('%s_marker_gene_dict_violin.png' % name)

sc.pl.dotplot(adata,marker_genes_dict,'leiden',dendrogram=True)
plt.tight_layout()

plt.savefig('%s_marker_gene_dict_dotplot.png' % name)

cluster2annotation0={
    #total 41: 0-40
    '27':'fibroblasts',
    '35':'fibroblasts',
    '25':'fibroblasts',
    '6':'fibroblasts',
    '31':'fibroblasts',
    '29':'fibroblasts',
    '0':'fibroblasts',
    '30':'fibroblasts',
    '10':'fibroblasts',
    '7':'fibroblasts',
    '16':'fibroblasts',




    '11':'Bcells',


    '21':'Macrophages', #check
    '17':'Macrophages', #check

    #'21':'myloid',


    '26':'cancercells',
    '4':'cancercells',
    '32':'cancercells', #important?
    '12':'cancercells',
    '3':'cancercells',
    '18':'cancercells',

    '19':'cancercells', 
    '14':'cancercells',
    '28':'cancercells',
    '13':'cancercells',
    '20':'cancercells',
    #'34':'cancercells', #check
    #'36':'myofibroblasts', #check





    '2':'Tcells',
    '23':'Tcells',
    '1':'Tcells',
    '5':'Tcells',




    '34':'plasma',
    '24':'plasma',
    '22':'plasma',


    '33':'endothelial',
    '15':'endothelial',
    '8':'endothelial',
    '9':'endothelial',





}


cluster2annotation=cluster2annotation0

for temp_item in range(len(cluster2annotation)):
    print(temp_item,cluster2annotation[str(temp_item)])

adata.obs['celltype']=adata.obs['leiden'].map(cluster2annotation).astype('category')
sc.pl.umap(adata,color='celltype',legend_loc='on data',frameon=False,legend_fontsize=12,legend_fontoutline=2,title='cell type',save='%s_leiden_celltype.png' % name)

adata.write('%s_population.h5ad' % name)

adata_new=adata[(adata.obs['celltype']=='CXCL14cancer')|(adata.obs['celltype']=='cancercells')|(adata.obs['celltype']=='cancerstemcells'),:]

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

