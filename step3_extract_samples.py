import scanpy as sc

name='celltype_name_bbknn_reduce_inhouse_preprocessed_nonorm.h5ad'

outname='celltype_name_bbknn_reduce_inhouse_preprocessed_nonorm_11samples.h5ad'

adata=sc.read(name)

remove_samples=['4','5','9']

new_adata=adata[~adata.obs['sample'].isin(remove_samples),:]
new_adata.write(outname)

#filter cell less than 200 expressed genes
sc.pp.filter_cells(new_adata,min_genes=200)
print(new_adata)