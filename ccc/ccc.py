# ccc
import scanpy as scc
import sys
import matplotlib.pyplot as plt
import os

name=sys.argv[1]
adata_concat=scc.read(name)

from liana.method import cellphonedb
import liana as li

# cellphone db

# for big type

cellphonedb(
    adata_concat, groupby="final_celltype", use_raw=False, return_all_lrs=True, verbose=True,layer="log1norm"
)

fig=li.pl.dotplot(
    adata=adata_concat,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    target_labels = ['bcell', 'endo', 'epithelial', 'fibroblast', 'macrophagecell',     'pericytes', 'plasma', 'tcell'],
    source_labels=["fibroblast"],
    # since cpdbv2 suggests using a filter to FPs
    # we can filter the interactions according to p-values <= 0.01
    filterby="cellphone_pvals",
    filter_lambda=lambda x: x <= 0.001,
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,  # we want to prioritize those with highest expression
    top_n=100,  # and we want to keep only the top 20 interactions
    figure_size=(15, 20),
    size_range=(1, 6),
)
fig.save('%s_fibroblast_epithelial_cellphone_celltype.png' % os.path.basename(name).split('.h5ad')[0])
adata_concat.uns["liana_res"].to_csv('%s_cellphone_celltype.csv' % os.path.basename(name).split('.h5ad')[0])

#for subtype
cellphonedb(
    adata_concat, groupby="final_subcelltype", use_raw=False, return_all_lrs=True, verbose=True,layer="log1norm"
)

fig=li.pl.dotplot(
    adata=adata_concat,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    target_labels =["Epithelial cells_0", "Epithelial cells_1", "Epithelial cells_2","Epithelial cells_3","Epithelial cells_4","Epithelial cells_5", "Epithelial cells_6", "Epithelial cells_7", "Epithelial cells_8"],
    source_labels=["Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2","Fibroblasts_3","Fibroblasts_4", "Fibroblasts_5","Fibroblasts_6","Fibroblasts_7" ],
    # since cpdbv2 suggests using a filter to FPs
    # we can filter the interactions according to p-values <= 0.01
    filterby="lr_means",
    filter_lambda=lambda x: x >1.5,
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    # expression magnitude
    orderby="cellphone_pvals",
    orderby_ascending=True,  # we want to prioritize those with highest expression
    top_n=100,  # and we want to keep only the top 20 interactions
    figure_size=(20, 25),
    size_range=(1, 6),
)
fig.save('%s_cellphone.png' % os.path.basename(name).split('.h5ad')[0])



fig=li.pl.dotplot(
    adata=adata_concat,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    target_labels =["T cells_0", "T cells_1", "T cells_2", "T cells_3", "T cells_4", "T cells_5",  "T cells_6",  "T cells_7",  "T cells_8", ],
    source_labels=["Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2","Fibroblasts_3","Fibroblasts_4", "Fibroblasts_5","Fibroblasts_6","Fibroblasts_7" ],
    # since cpdbv2 suggests using a filter to FPs
    # we can filter the interactions according to p-values <= 0.01
    filterby="cellphone_pvals",
    filter_lambda=lambda x: x <= 0.01,
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    # expression magnitude
    orderby="cellphone_pvals",
    orderby_ascending=True,  # we want to prioritize those with highest expression
    top_n=100,  # and we want to keep only the top 20 interactions
    figure_size=(20, 15),
    size_range=(1, 6),
)
fig.save('%s_fibroblast_tcells_cellphone.png' % os.path.basename(name).split('.h5ad')[0])


# B cells
b_cells_label=set(adata_concat.obs['final_subcelltype'][adata_concat.obs['final_subcelltype'].str.contains('B cells')])

fig=li.pl.dotplot(
    adata=adata_concat,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    target_labels =list(b_cells_label),
    source_labels=["Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2","Fibroblasts_3","Fibroblasts_4", "Fibroblasts_5","Fibroblasts_6","Fibroblasts_7" ],
    # since cpdbv2 suggests using a filter to FPs
    # we can filter the interactions according to p-values <= 0.01
    filterby="cellphone_pvals",
    filter_lambda=lambda x: x <= 0.01,
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    # expression magnitude
    orderby="cellphone_pvals",
    orderby_ascending=True,  # we want to prioritize those with highest expression
    top_n=100,  # and we want to keep only the top 20 interactions
    figure_size=(20, 25),
    size_range=(1, 6),
)
fig.save('%s_fibroblast_bcells_cellphone.png' % os.path.basename(name).split('.h5ad')[0])

#macrophage
macrophage_label=set(adata_concat.obs['final_subcelltype'][adata_concat.obs['final_subcelltype'].str.contains('Macrophages')])
fig=li.pl.dotplot(
    adata=adata_concat,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    target_labels =list(macrophage_label),
    source_labels=["Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2","Fibroblasts_3","Fibroblasts_4", "Fibroblasts_5","Fibroblasts_6","Fibroblasts_7" ],
    # since cpdbv2 suggests using a filter to FPs
    # we can filter the interactions according to p-values <= 0.01
    filterby="cellphone_pvals",
    filter_lambda=lambda x: x <= 0.01,
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,  # we want to prioritize those with highest expression
    top_n=100,  # and we want to keep only the top 20 interactions
    figure_size=(20, 25),
    size_range=(1, 6),
)
fig.save('%s_fibroblast_macrophage_cellphone.png' % os.path.basename(name).split('.h5ad')[0])


adata_concat.uns["liana_res"].to_csv('%s_cellphone_subcelltype.csv' % os.path.basename(name).split('.h5ad')[0])


fig=li.pl.dotplot(
    adata=adata_concat,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    source_labels=["T cells_0", "T cells_1", "T cells_2", "T cells_3", "T cells_4", "T cells_5",  "T cells_6",  "T cells_7",  "T cells_8", ],
    target_labels=["Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2","Fibroblasts_3","Fibroblasts_4", "Fibroblasts_5","Fibroblasts_6","Fibroblasts_7" ],
    # since cpdbv2 suggests using a filter to FPs
    # we can filter the interactions according to p-values <= 0.01
    filterby="cellphone_pvals",
    filter_lambda=lambda x: x <= 0.01,
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    # expression magnitude
    orderby="cellphone_pvals",
    orderby_ascending=True,  # we want to prioritize those with highest expression
    top_n=100,  # and we want to keep only the top 20 interactions
    figure_size=(20, 15),
    size_range=(1, 6),
)
fig.save('%s_tcells_fibroblast_cellphone.png' % os.path.basename(name).split('.h5ad')[0])

fig=li.pl.dotplot(
    adata=adata_concat,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot

    #exchange target_labels and source_labels
    target_labels =["Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2","Fibroblasts_3","Fibroblasts_4", "Fibroblasts_5","Fibroblasts_6","Fibroblasts_7" ],
    source_labels =["Epithelial cells_0", "Epithelial cells_1", "Epithelial cells_2","Epithelial cells_3","Epithelial cells_4","Epithelial cells_5", "Epithelial cells_6", "Epithelial cells_7", "Epithelial cells_8"],
    # since cpdbv2 suggests using a filter to FPs
    # we can filter the interactions according to p-values <= 0.01
    filterby="cellphone_pvals",
    filter_lambda=lambda x: x <= 0.01,
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    # expression magnitude
    orderby="cellphone_pvals",
    orderby_ascending=True,  # we want to prioritize those with highest expression
    top_n=100,  # and we want to keep only the top 20 interactions
    figure_size=(20, 15),
    size_range=(1, 6),
)

fig.save('%s_cellphone2.png' % os.path.basename(name).split('.h5ad')[0])

from liana.method import rank_aggregate

rank_aggregate(
    adata_concat, groupby="final_subcelltype", return_all_lrs=True, use_raw=False, verbose=True
)

adata_concat.uns["liana_res"].drop_duplicates(
    ["ligand_complex", "receptor_complex"]
).head()

plt.close()
fig=li.pl.dotplot(
    adata=adata_concat,
    colour="magnitude_rank",
    size="specificity_rank",
    inverse_colour=True,  # we inverse sign since we want small p-values to have large sizes
    inverse_size=True,
    # We choose only the cell types which we wish to plot
    target_labels =["Epithelial cells_0", "Epithelial cells_1", "Epithelial cells_2","Epithelial cells_3","Epithelial cells_4","Epithelial cells_5", "Epithelial cells_6", "Epithelial cells_7", "Epithelial cells_8"],
    source_labels=["Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2","Fibroblasts_3","Fibroblasts_4", "Fibroblasts_5","Fibroblasts_6","Fibroblasts_7" ],
 
    # since the rank_aggregate can also be interpreted as a probability distribution
    # we can again filter them according to their specificity significance
    # yet here the interactions are filtered according to
    # how consistently highly-ranked is their specificity across the methods
    filterby="specificity_rank",
    filter_lambda=lambda x: x <= 0.05,
    # again, we can also further order according to magnitude
    orderby="magnitude_rank",
    orderby_ascending=True,  # prioritize those with lowest values
    top_n=30,  # and we want to keep only the top 20 interactions
    figure_size=(15, 20),
    size_range=(1, 6),
    return_fig=True,
)

#ggplot save
fig.save('%s_liana.png' % os.path.basename(name).split('.h5ad')[0],dpi=100)