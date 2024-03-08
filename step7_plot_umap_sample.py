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
import matplotlib.font_manager as fm

#font arial: /disk3/users/lvwen/ARIAL.TTF
#local font file
font_path='/disk3/users/lvwen/ARIAL.TTF'
myfont=fm.FontProperties(fname=font_path)
name=sys.argv[1] # input name
#gene_name=sys.argv[2]

adata=sc.read(name)

#marker_genes=['DCN','GSN','APOD','CXCL12','CLDN4','TACSTD2','TM4SF1','MGP','NKG7','B2M','HLA-B','HLA-A','VIM', 'CRYAB','S100B','LGALS1']
#marker_genes=open(gene_name).read().split()


#vp=sc.pl.stacked_violin(adata,marker_genes,groupby='clusters',rotation=90,return_fig=True,dendrogram=True)
#vp=sc.pl.umap(adata,color=new_genes,return_fig=True)
#add string sample before sample name
adata.obs['sample']='sample_'+adata.obs['sample'].astype(str)

#14 samples, 7 PT, sample_0: indiviudal1_T; sample_1: PT1_N; sample_2: PT2_T; sample_3: PT2_N; sample_4: PT3_T; sample_5: PT3_N; sample_6: PT4_T; sample_7: PT4_N; sample_8: PT5_T; sample_9: PT5_N; sample_10: PT6_T; sample_11: PT6_N; sample_12: PT7_T; sample_13: PT7_N
adata.obs['sample_new']=adata.obs['sample'].replace({'sample_0':'PT1_T','sample_1':'PT1_N','sample_2':'PT2_T','sample_3':'PT2_N','sample_4':'PT3_T','sample_5':'PT3_N','sample_6':'PT4_T','sample_7':'PT4_N','sample_8':'PT5_T','sample_9':'PT5_N','sample_10':'PT6_T','sample_11':'PT6_N','sample_12':'PT7_T','sample_13':'PT7_N'})
vp=sc.pl.umap(adata,color='sample_new',return_fig=True,frameon=False,legend_fontsize=5,size=0.5 ,title=None)
vp.set_size_inches(2.5,2)
#get fig title

plt.title(None,fontproperties=myfont,fontsize=5,fontweight='bold')
plt.xlabel('UMAP1',fontproperties=myfont,fontsize=5)
plt.ylabel('UMAP2',fontproperties=myfont,fontsize=5)
#set title
#plt.title('Samples',fontproperties=myfont,fontsize=5)

#plt.savefig('%s_umap.png' % sample)
#vp.add_totals().style(ylim=(0,8))
plt.tight_layout()

vp.savefig('step7_%s_%s_marker_genes_umap.png' % ('sample',os.path.basename(name)),dpi=300)
vp.clf()



#stacked bar chart for cell types percentage in each sample
# figure size 6,6
adata.obs['sample']=adata.obs['sample_new']
samples=adata.obs['sample'].unique()
celltypes=adata.obs['celltype_new'].unique()
celltype_percent=pd.DataFrame(index=samples,columns=celltypes)
for sample in samples:
    for celltype in celltypes:
        celltype_percent.loc[sample,celltype]=adata.obs[(adata.obs['sample']==sample) & (adata.obs['celltype_new']==celltype)].shape[0]
celltype_percent=celltype_percent.fillna(0)
celltype_percent=celltype_percent.astype(int)
celltype_percent=celltype_percent.div(celltype_percent.sum(axis=1),axis=0)
celltype_percent=celltype_percent.fillna(0)
celltype_percent=celltype_percent.round(3)
celltype_percent=celltype_percent.sort_index(axis=1)
celltype_percent=celltype_percent.sort_index(axis=0)
celltype_percent.to_csv('step7_%s_%s_celltype_percent.csv' % ('sample',os.path.basename(name)))
#plot stacked bar chart, for each sample
fig = plt.figure()

ax=celltype_percent.plot.bar(stacked=True,figsize=(2.5,2),fontsize=5,legend=True)
ax.grid(False)
ax.legend(fontsize=5,loc='center left', bbox_to_anchor=(1, 0.5),frameon=False)
#change legend size
#plt.legend(loc='upper right',prop=myfont,fontsize="1")
plt.xticks(rotation=90,fontproperties=myfont,fontsize=5)
plt.yticks(fontproperties=myfont,fontsize=5)
plt.ylabel('Percentage',fontproperties=myfont,fontsize=5)
#plt.xlabel('Samples',fontproperties=myfont,fontsize=5)
plt.title('Cell type percentage',fontproperties=myfont,fontsize=5)
plt.tight_layout()
plt.savefig('step7_%s_%s_celltype_percent.png' % ('sample',os.path.basename(name)),dpi=300)
plt.clf()


#change celltype name according dict
print(adata.obs['celltype_new'].unique())
names_new={'Tcells':'T cells','myofibroblasts':'Fibroblasts','cancercells':'Epithelial cells','endothelial':'Endothelial cells','pericytes':'Pericytes','plasma':'Plasma cells','Macrophages':'Macrophages','Bcells':'B cells'}

adata.obs['celltype_new']=adata.obs['celltype_new'].replace(names_new)
vp=sc.pl.umap(adata,color='celltype_new',add_outline=False,legend_loc='on data',use_raw=False,legend_fontoutline=1,frameon=True,legend_fontsize=5,size=2,return_fig=True, )
vp.set_size_inches(2,2)

#title font size
plt.title('Cell type',fontproperties=myfont,fontsize=5,fontweight='bold')
# xlabel size
plt.xlabel('UMAP1',fontproperties=myfont,fontsize=5)
plt.ylabel('UMAP2',fontproperties=myfont,fontsize=5)

#plt.savefig('%s_umap.png' % sample)
#vp.add_totals().style(ylim=(0,8))
plt.tight_layout()

vp.savefig('step7_%s_%s_marker_genes_umap.png' % ('celltype',os.path.basename(name)),dpi=300)
plt.close()


vp=sc.pl.umap(adata,color='leiden_res0_4',add_outline=False,legend_loc='on data',use_raw=False,legend_fontoutline=1,frameon=True,legend_fontsize=5,size=2,return_fig=True, )
vp.set_size_inches(2,2)

#title font size
plt.title('Leiden',fontproperties=myfont,fontsize=5,fontweight='bold')
# xlabel size
plt.xlabel('UMAP1',fontproperties=myfont,fontsize=5)
plt.ylabel('UMAP2',fontproperties=myfont,fontsize=5)

#plt.savefig('%s_umap.png' % sample)
#vp.add_totals().style(ylim=(0,8))
plt.tight_layout()

vp.savefig('step7_%s_%s_marker_genes_umap.png' % ('leiden',os.path.basename(name)),dpi=300)
plt.close()


#umap for marker genes

#marker_genes=['DCN','GSN','APOD','CXCL12','CLDN4','TACSTD2','TM4SF1','MGP','NKG7','B2M','HLA-B','HLA-A','VIM', 'CRYAB','S100B','LGALS1']
import matplotlib as mpl
gene_name='marker_genes'
marker_genes=open('marker_genes').read().split()

new_genes=[]
for gene in marker_genes:
    if gene in adata.var.index:
        new_genes.append(gene)

print(set(marker_genes)-set(new_genes))
#vp=sc.pl.stacked_violin(adata,marker_genes,groupby='clusters',rotation=90,return_fig=True,dendrogram=True)
#vp=sc.pl.umap(adata,color=new_genes,return_fig=True)

# cmap gray to red
cmap=mpl.colors.LinearSegmentedColormap.from_list('mycmap', ['gray', 'red'])
vp=sc.pl.umap(adata,color=marker_genes,return_fig=True,ncols=3,cmap=cmap,vmin=0, add_outline=True)

#plt.savefig('%s_umap.png' % sample)
#vp.add_totals().style(ylim=(0,8))
vp.savefig('step7_%s_%s_marker_genes_umap.png' % (gene_name,os.path.basename(name)),dpi=300)