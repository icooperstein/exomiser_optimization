import pandas as pd
import altair as alt
import os


patient_ids = [] ##list of patient IDs Exomiser (or Genomiser) was run on

all_results_df = pd.DataFrame()
path = 'your/path/Exomiser/exomiser_results/'
RUN_TYPES=['exomiser_default', 'filtered_exomiser_default']

for ID in patient_ids:
    for RUN_TYPE in RUN_TYPES:
        output_filename = str(ID)+'/results/' +str(ID) + "_" + RUN_TYPE+ '.genes.tsv'
        if os.path.exists(path+output_filename):
            gene_result = pd.read_csv(path+output_filename, sep='\t')
            print(len(gene_result), 'genes')
            gene_result['Patient_ID'] = [ID]*len(gene_result)

            gene_result['RunType'] = [str(RUN_TYPE)] * len(gene_result)
            patient_df = gene_result[['Patient_ID', '#RANK', 'GENE_SYMBOL', 'P-VALUE', 'EXOMISER_GENE_COMBINED_SCORE', 'EXOMISER_GENE_PHENO_SCORE', 'EXOMISER_GENE_VARIANT_SCORE', 'RunType']]
            running_results = pd.concat([all_results_df, patient_df])
            all_results_df = running_results.copy()
        else:
            print(path+output_filename)


all_results_df=all_results_df.reset_index(drop=True)


###save results table
all_results_df.to_csv('complete_exomiser_results.tsv', sep='\t', index=None)