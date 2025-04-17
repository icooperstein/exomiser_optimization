import pandas as pd
import altair as alt
import os


patient_ids = [] ##list of patient IDs Exomiser (or Genomiser) was run on

all_results_df = pd.DataFrame()
path = 'your/path/Exomiser/exomiser_results/'
RUN_TYPES=['exomiser_default', 'filtered_exomiser_default']

for ID in patient_ids:
    #if file in casesToCheck:
    for RUN_TYPE in RUN_TYPES:
        variant_file =  str(ID)+'/2406_results_newCADD/' +str(ID) + "_" + RUN_TYPE+ '.variants.tsv'
        if os.path.exists(path+variant_file):
            variantResult = pd.read_csv(path+variant_file, sep='\t')
            variantResult = variantResult[variantResult['CONTRIBUTING_VARIANT']==1]
            print(len(variantResult), 'contributing variants')
            variantResult['Patient_ID'] = [ID]*len(variantResult)
            variantResult['RunType'] = [str(RUN_TYPE)] * len(variantResult)
            patient_df = variantResult[['Patient_ID', '#RANK', 'ID','GENE_SYMBOL', 'P-VALUE', 'EXOMISER_GENE_COMBINED_SCORE', 'EXOMISER_GENE_PHENO_SCORE', 'EXOMISER_GENE_VARIANT_SCORE', 'EXOMISER_VARIANT_SCORE', 'MAX_PATH_SOURCE','MAX_PATH', 'ALL_PATH','RunType']]
            running_results = pd.concat([all_results_df, patient_df])
            all_results_df = running_results.copy()
        else:
            print(path+output_filename)


all_results_df=all_results_df.reset_index(drop=True)

###save results table
all_results_df.to_csv('complete_exomiser_results.tsv', sep='\t', index=None)