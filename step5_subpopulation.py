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



def ananlysis_subset_data(adatasub,sample):
    #adata=adata_merge
    #adata.write('step3_res_%s_subpopulation_raw.h5ad' % sample)
    adatasub.X=adatasub.layers['counts'].toarray()
    
    print(adatasub)
    sc.pp.filter_genes(adatasub,min_cells=3)
    print(adatasub)
    sc.pp.normalize_total(adatasub,target_sum=1e4)
    sc.pp.log1p(adatasub)

    resolution_value=0.8 #0.5 OK
    sc.pp.pca(adatasub)
    sc.pp.neighbors(adatasub)
    sc.tl.umap(adatasub)
    sc.tl.leiden(adatasub, resolution=resolution_value)

    sc.pl.umap(adatasub, color=['sample','leiden'])

    plt.savefig('%s_test1.png' % sample,dpi=100)
    sc.tl.rank_genes_groups(adatasub,'leiden',methods='t-test' ,pts=True)
    sc.pl.rank_genes_groups(adatasub,n_genes=25,sharey=False)
    plt.tight_layout()

    plt.savefig('%s_rank_genes_t-test.png' % sample,dpi=200)

    sc.tl.rank_genes_groups(adatasub,'leiden',methods='wilcoxon' ,pts=True)
    sc.pl.rank_genes_groups(adatasub,n_genes=25,sharey=False)
    plt.tight_layout()

    plt.savefig('%s_rank_genes_wilcoxon.png' % sample,dpi=200)

    bbknn.bbknn(adatasub, batch_key='sample',neighbors_within_batch=2)
    sc.tl.umap(adatasub)
    sc.tl.leiden(adatasub, resolution=resolution_value)

    sc.pl.umap(adatasub, color=['sample','leiden'])
    plt.savefig('%s_test4.png' % sample ,dpi=100)
    sc.tl.rank_genes_groups(adatasub,'leiden',methods='t-test' ,pts=True)
    sc.pl.rank_genes_groups(adatasub,n_genes=25,sharey=False)
    plt.tight_layout()

    plt.savefig('%s_rank_genes_bbknn_t-test.png' % sample,dpi=200)

    sc.tl.rank_genes_groups(adatasub,'leiden',methods='wilcoxon' ,pts=True)
    sc.pl.rank_genes_groups(adatasub,n_genes=25,sharey=False)
    plt.tight_layout()

    plt.savefig('%s_rank_genes_bbknn_wilcoxon.png' % sample,dpi=200)

    if not isinstance(adatasub.X, np.ndarray):
        adatasub.X=adatasub.X.toarray()
    bbknn.ridge_regression(adatasub, batch_key=['sample'], confounder_key=['leiden'])
    sc.pp.pca(adatasub)
    bbknn.bbknn(adatasub, batch_key='sample',neighbors_within_batch=2)
    sc.tl.umap(adatasub)
    sc.tl.leiden(adatasub, resolution=resolution_value)

    sc.pl.umap(adatasub, color=['sample','leiden'])
    plt.savefig('%s_test5.png' % sample,dpi=100)

    #adata.var['mt']=adata.var_names.str.startswith('MT-')
    #sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,log1p=False,inplace=True)
    #sc.pp.highly_variable_genes(adata,min_mean=0.0125,max_mean=3,min_disp=0.5)
    #sc.pp.highly_variable_genes(adata,min_disp=0.5)

    #sample='merged_data' #setting the sample id, 
    #sc.pl.highly_variable_genes(adata)
    #plt.tight_layout()
    #plt.savefig('%s_highly_variable_genes.png' % sample)

    #adata.raw=adata
    #adata=adata[:,adata.var.highly_variable]
    #print(adata)
    #sc.pp.regress_out(adata,['total_counts','pct_counts_mt']) # need or not?
    #sc.pp.scale(adata,max_value=10)
    #sc.pp.regress_out(adata,['total_counts','pct_counts_mt']) # need or not?
    #sc.pp.scale(adata,max_value=10)

    #sc.tl.pca(adata,svd_solver='arpack')
    #sc.pl.pca(adata,color='SFRP4')
    #plt.savefig('%s_SFRP4_PCA.png' % sample)

    sc.pl.pca_variance_ratio(adatasub,log=True)
    plt.savefig('%s_variance_ratio.png' % sample)


    #sc.pp.neighbors(adata,n_neighbors=20,n_pcs=40)
    #sc.tl.umap(adata)
    sc.pl.umap(adatasub,color=marker_genes)
    plt.tight_layout()

    plt.savefig('%s_umap.png' % sample)

    #sc.tl.leiden(adata,key_added='clusters',resolution=0.3)
    sc.pl.umap(adatasub,color='leiden',add_outline=False,legend_loc='on data',use_raw=False, legend_fontsize=12,legend_fontoutline=2,frameon=True,title='clustering of cells',palette='Set1')
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.axis()
    plt.tight_layout()

    plt.savefig('%s_umap_leiden.png' % sample ,dpi=200)

    sc.tl.rank_genes_groups(adatasub,'leiden',methods='t-test' ,pts=True,layer='log1norm')
    sc.pl.rank_genes_groups(adatasub,n_genes=25,sharey=False)
    plt.tight_layout()

    plt.savefig('%s_rank_genes_t-test.png' % sample,dpi=200)

    #sc.tl.rank_genes_groups(adata,'clusters',methods='t-test',pts=True)
    #sc.pl.rank_genes_groups(adata,n_genes=25,sharey=False)
    #plt.savefig('%s_rank_genes_ttest.png' % sample,dpi=200)

    vp=sc.pl.stacked_violin(adatasub,marker_genes,groupby='leiden',rotation=90,return_fig=True)
    vp.add_totals().style(ylim=(0,8))
    vp.savefig('%s_marker_genes_volin.png' % sample,dpi=300)

    sc.tl.dendrogram(adatasub,groupby='leiden')
    sc.pl.rank_genes_groups_dotplot(adatasub,n_genes=6)
    #plt.tight_layout()

    plt.savefig('%s_dotplot.png' % sample,dpi=200)

    #sc.pl.rank_genes_groups_heatmap(adata,n_genes=10,use_raw=False, swap_axes=True,vmin=-3,vmax=3,cmap='bwr',figsize=(15,5),show=False)
    #plt.savefig('%s_heatmap.png' % sample,dpi=200)

    sc.pl.rank_genes_groups_stacked_violin(adatasub,n_genes=6,cmap='viridis_r')
    #plt.tight_layout()

    plt.savefig('%s_stacked_violin.png' % sample,dpi=200)

    #sc.pl.tracksplot(adata,marker_genes_dict,groupby='clusters',dendrogram=True)
    #plt.savefig('%s_tracksplot.png' % sample,dpi=200)

    adatasub.write('step3_res_%s_subpopulation.h5ad' % sample)

    result=adatasub.uns['rank_genes_groups']
    groups=result['names'].dtype.names
    temp_data=pd.DataFrame({group+'_'+key: result[key][group] for group in groups for key in ['names','pvals','pvals_adj','logfoldchanges','pts','pts_rest']})
    temp_data.to_csv('step3_%s_subpopulation_markers.csv'  % sample)




marker_genes_dict={'NK':['GNLY','NKG7'],
'T-cell':['CD3D'],
'B-cell':['CD79A','MS4A1'],
'Fibroblast':['DCN'],}
#'epithelial':['KRT18','XBP1']}

marker_genes_dict={
    'cancercells':['KRT19'],
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

samples=['X1Ca','X1Pa','X2Ca','X2Pa','X3Ca','X3Pa','X4Ca','X4Pa','X5Ca','X5Pa','X6Ca','X6Pa','X7Ca','X7Pa']

marker_genes=['DCN','CFD','COL1A2','MMP2']

import sys

name=sys.argv[1]
name='celltype_name_bbknn_reduce_inhouse_preprocessed_nonorm_11samples.h5ad_population.h5ad'
adata=sc.read(name)

#if 'log1p' in adata_merge.uns:
    #adata_merge.X=np.expm1(adata_merge.X)
#    adata_merge.uns['log1p']["base"] = None
if 'base' not in adata.uns['log1p']:
    adata.uns['log1p']["base"] = None

#adata_fibroblast=adata[adata.obs['celltype_new'].isin(['Fibroblasts','Perivascular cells'])] # fibroblast
adata_fibroblast=adata[adata.obs['celltype_new'].isin(['Fibroblasts'])] # fibroblast

adata_epi=adata[adata.obs['celltype_new'].isin(['Epithelial cells'])]
adata_tcell=adata[adata.obs['celltype_new'].isin(['T cells'])]
#adata_bcell=adata[adata.obs['clusters'].isin(['25','12','10','28'])]
adata_bcell=adata[adata.obs['celltype_new'].isin(['B cells'])]
adata_pericytes=adata[adata.obs['celltype_new'].isin(['Perivascular cells'])]
adata_macrophagecell=adata[adata.obs['celltype_new'].isin(['Macrophages'])]
adata_endo=adata[adata.obs['celltype_new'].isin(['Endothelial cells'])]
adata_plasma=adata[adata.obs['celltype_new'].isin(['Plasma cells'])]
#adata_myeloid=adata[adata.obs['celltype'].isin(['myloid'])] #wrong for myeloid
#adata_tme=adata[adata.obs['celltype_new'].isin(['plasma','fibroblasts'])]

#adata_new=adata[(adata.obs['celltype']=='CXCL14cancer')|(adata.obs['celltype']=='cancercells')|(adata.obs['celltype']=='cancerstemcells'),:]




ananlysis_subset_data(adata_fibroblast,'fibroblast_bbknn_regression')
ananlysis_subset_data(adata_epi,'epithelial_bbknn_regression')
ananlysis_subset_data(adata_tcell,'tcell_bbknn_regression')
ananlysis_subset_data(adata_bcell,'bcell_bbknn_regression')
ananlysis_subset_data(adata_pericytes,'pericytes_bbknn_regression')
ananlysis_subset_data(adata_macrophagecell,'macrophagecell_bbknn_regression')
ananlysis_subset_data(adata_endo,'endo_bbknn_regression')
ananlysis_subset_data(adata_plasma,'plasma_bbknn_regression')
#ananlysis_subset_data(myeloid,'myeloid_bbknn_regression')
