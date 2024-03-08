# usage: python3 /home/lu/pipeline/tools/gseapy/run_over_representation_analysis_sc_csv.py /home/lu/projects/single_cell_breast_cancer/final/step3_fibroblast_bbknn_regression_subpopulation_markers.csv
# run all
import sys
import os
import pandas as pd

csv_name = sys.argv[1]
data=pd.read_csv(csv_name,sep=',',index_col=0)

top_gene_num=1500
subtypes=list(data.columns[data.columns.str.contains('_names')])
for subtype in subtypes:
    #if subtype!='1_names':
    #    print(subtype)
    #    continue
    basename=os.path.basename(csv_name).split('.')[0]
    prefix=subtype.split('_')[0]

    # Path: tools/gseapy/run_over_representation_analysis_batch.py

    enrich_folder=basename+'_%s_enrichr' % subtype

    if not os.path.exists(enrich_folder):
        os.makedirs(enrich_folder)

    gene_name=enrich_folder+'_genes'

    sig_set=data[data['%s_pvals_adj' % prefix]<0.05]
    sig_set_up=sig_set[sig_set['%s_logfoldchanges' % prefix]>0]
    sig_set_down=sig_set[sig_set['%s_logfoldchanges' % prefix]<0]
    #sort by column logfoldchanges
    sig_set_up.sort_values(ascending=False,inplace=True,by='%s_logfoldchanges' % prefix)
    sig_set_down.sort_values(ascending=True,inplace=True,by='%s_logfoldchanges' % prefix)

    out=open(gene_name+'_up','w')
    for item in sig_set_up.index[:top_gene_num]:
        out.write(sig_set_up.loc[item][subtype]+'\n')
    out.close()
    out=open(gene_name+'_down','w')
    for item in sig_set_down.index[:top_gene_num]:
        out.write(sig_set_down.loc[item][subtype]+'\n')
    out.close()

    #for prerank genes according to logfoldchanges
    sig_set.sort_values(ascending=False,inplace=True,by='%s_logfoldchanges' % prefix)
    out=open(gene_name+'_prerank','w')
    for item in sig_set.index:
        out.write(sig_set.loc[item][subtype]+'\t'+str(sig_set.loc[item]['%s_logfoldchanges' % prefix])+'\n')
    out.close()

    gene_name_up=gene_name+'_up'
    gene_name_down=gene_name+'_down'
    #check

    #print('gseapy enrichr -i %s -g GO_Molecular_Function_2021 -v -o %s/enrichr_MF_%s' %(gene_name, enrich_folder, basename))
    basename=basename+'_up'
    enrich_folder=basename+'_%s_enrichr_up' % subtype

    if not os.path.exists(enrich_folder):
        os.makedirs(enrich_folder)

    os.system('gseapy enrichr -i %s -g GO_Molecular_Function_2021 -v -o %s/enrichr_MF_%s' %(gene_name_up, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g GO_Cellular_Component_2021 -v -o %s/enrichr_CC_%s' %(gene_name_up, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g GO_Biological_Process_2021 -v -o %s/enrichr_BP_%s' %(gene_name_up, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g KEGG_2021_Human -v -o %s/enrichr_KEGG_%s' %(gene_name_up, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g Reactome_2022 -v -o %s/enrichr_Reactome_%s' %(gene_name_up, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g miRTarBase_2017 -v -o %s/enrichr_miRTarBase_%s' %(gene_name_up, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g MSigDB_Oncogenic_Signatures -v -o %s/enrichr_MSigDB_Oncogenic_%s' %(gene_name_up, enrich_folder, basename))

    basename=basename+'_down'
    enrich_folder=basename+'_%s_enrichr_down' % subtype

    if not os.path.exists(enrich_folder):
        os.makedirs(enrich_folder)

    os.system('gseapy enrichr -i %s -g GO_Molecular_Function_2021 -v -o %s/enrichr_MF_%s' %(gene_name_down, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g GO_Cellular_Component_2021 -v -o %s/enrichr_CC_%s' %(gene_name_down, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g GO_Biological_Process_2021 -v -o %s/enrichr_BP_%s' %(gene_name_down, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g KEGG_2021_Human -v -o %s/enrichr_KEGG_%s' %(gene_name_down, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g Reactome_2022 -v -o %s/enrichr_Reactome_%s' %(gene_name_down, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g miRTarBase_2017 -v -o %s/enrichr_miRTarBase_%s' %(gene_name_down, enrich_folder, basename))
    os.system('gseapy enrichr -i %s -g MSigDB_Oncogenic_Signatures -v -o %s/enrichr_MSigDB_Oncogenic_%s' %(gene_name_down, enrich_folder, basename))


