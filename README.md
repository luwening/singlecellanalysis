# singlecellanalysis
The code is used for reproduce for publicatioin: XXXX

**scRNA-seq data analysis**

Generally, it can divided into following steps:
1. merge 14 samples into a matrix file

python step1_merge_data.py

2. preprocess the matrix

python step2_preprocess.py

3. extract samples (remove samples with low quality) and bbknn

python step3_extract_samples.py

python step3_merge_bbknn.py #bbknn

4. cell population clustering (cell type level)

python step4_populations.py

5. CAF subpopulation clustering (sub cell type level)

python step5_subpopulation.py

6. scripts used for visualization the single cell data

python step6_plot_genes_umap.py

python step6_plot_umap_sample.py

python step6_plot_genes_dotplot_conditions.py

**paga analysis**

python do_paga.py

**spatial analysis part**

python spatial.py

**enrichment analysis**

python run_over_representation_analysis_sc_csv.py

python prerank_test.py

**conda enviroment**

details in packages_in_conda.txt

