# singlecellanalysis
The code is used for reproduce for publicatioin: XXXX

inhouse scRNA-seq data analysis

Generally, it can divided into following steps:
1. merge 14 samples into a matrix file

python step1_merge_data.py

2. preprocess the matrix

python step2_preprocess.py

3. extract samples (remove samples with low quality)

python step3_extract_samples.py

2. cell population clustering (cell type level)

3. CAF subpopulation clustering (sub cell type level)

4. scripts used for visualization the single cell data


spatial analysis part
python spatial.py

#conda enviroment

