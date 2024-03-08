# working path: /disk3/users/lvwen/download/10x/
# usage: python3 ~/pipeline/single_cell/10x/scanpy_pipeline/spatial/test.py

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

adata=sc.read_visium(path="/disk3/users/lvwen/download/10x/",count_file='V1_Breast_Cancer_Block_A_Section_1_raw_feature_bc_matrix.h5',library_id='V1_Breast_Cancer_Block_A_Section_1')
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
print(adata)
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
plt.savefig("test.png")

# Filter cells
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20]
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
plt.savefig("test2.png")

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
plt.savefig("test3.png")

sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)
plt.savefig("test4.png")

sc.pl.spatial(adata, img_key="hires", color="clusters", groups=["0", "5"], crop_coord=[1200, 1700, 1900, 1000], alpha=0.5, size=1.3)

sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="5", n_genes=10, groupby="clusters")


sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2","SFRP4",'SFRP1','SFRP2','WNT3A',"CFD",'KRT14','EPCAM','DCN','CD79A','CD3D','ERBB2'])
plt.savefig("test5.png")

sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2","SFRP4",'SFRP1','SFRP2','WNT3A',"CFD",'KRT14','EPCAM','DCN','CD79A','CD3D','ERBB2'], \
               groups = ["0"],crop_coord = [10000,20000,19000,10000], alpha = .9, size = 1.5)
plt.savefig("test5.png")

groups = ["4","8"],

#wnt pathway
sc.pl.spatial(adata, img_key="hires", color=["clusters", "WNT2",'WNT3','WNT6','WNT5A'])
plt.savefig("test6.png")

#wnt path
sc.pl.spatial(adata, img_key="hires", color=["clusters", "WNT2",'WNT3','WNT6','WNT5A','WNT3A'])
plt.savefig("test7.png")

#maker genes
sc.pl.spatial(adata, img_key="hires", color=["clusters", "WNT2",'WNT3','WNT6','WNT5A','WNT3A'])
plt.savefig("test7.png")

sc.tl.rank_genes_groups(adata, "clusters", inplace = True)
sc.pl.rank_genes_groups_heatmap(adata, groups = "4", groupby = "clusters", show = False)
output_name_prefix='spatial_test'
rank_genes_groups=adata.uns['rank_genes_groups']
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

gene_markers=[]
for cell_type in marker_genes_dict:
    gene_markers.extend(marker_genes_dict[cell_type])

#maker genes
sc.pl.spatial(adata, img_key="hires", color=["clusters", *gene_markers])
plt.savefig("test8.png")

sc.pl.dotplot(adata,marker_genes_dict,'clusters',dendrogram=True,figsize=(16,8))
plt.tight_layout(pad=10,h_pad=8)
#fig.subplots_adjust(left=0.2, right=0.98, top=0.9, bottom=0.1)
plt.savefig('%s_marker_gene_dict_dotplot.png' % output_name_prefix)