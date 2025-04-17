import pandas as pd
import altair as alt
'''
BENCHMARKING STRATEGY: Gene-level success (#1)
This is the least stringent success criteria.
The diagnostic gene is prioritized in the Exomiser/Genomiser output, regardless of the contributing variant's position or nucleotide change

ie: 
chr1:123456 A>T in Gene A
MOI: AR
exomiser form: 1-123456-A-T_AR

exomiser result: 1-56558-A-T_AR in Gene A

This is a SUCCESS by this metric because the diagnostic gene is in the output, regardless of what the variant is.

Necessary columns in input table:
1. ID (form patientID_gene)
2. genomic_coordinates_hg38 (chr1:123456 A>T)
3. MOI (AR, AD, XR, or XD)
Note: if performing singleton analysis, all Exomiser results will be in the form 1-123456-A-T_ANY so all variants will fail unless you set MOI column to ANY 

'''

variant_table = pd.read_csv('input.txt', sep='\t')
print(len(set(variant_table['ID'])),'genes')
print(len(variant_table), 'variants')



'''Define run-types'''

##FIGURE 1
RUN_TYPES = ['exomiser_default','filtered_exomiser_default']

'''Do the main thing'''
from collections import defaultdict
import math
data = defaultdict(list)
for column in variant_table.columns:
    data[column] = list(variant_table[column])[:]
failed=[]
for RUN_TYPE in set(RUN_TYPES):
    print(RUN_TYPE)
    n=0

    for i, row in variant_table.iterrows():
        
        full_ID = row['ID']
        gene = full_ID.split('_')[1]
        ID = full_ID.split('_')[0]

        try:
            results_file = pd.read_csv('/path/to/Exomiser/directory/exomiser_results/' + str(ID) +'/results/' + str(ID) + '_' + str(RUN_TYPE)+'.variants.tsv', sep='\t')
        except:
            print('no results file found for', str(ID))
            failed.append(ID)

        results_file = results_file[results_file['CONTRIBUTING_VARIANT']==1]

        parse_str = lambda x: x.split('_')[0]

        if gene in list(results_file['GENE_SYMBOL']):
            n+=1
            data['Gene_Level_' +str(RUN_TYPE)].append('Success')
            gene_result_table = results_file[results_file['GENE_SYMBOL']==gene]
            gene_result_rank = list(gene_result_table['#RANK'])[0]
            data['Gene_Level_rank_' +str(RUN_TYPE)].append(gene_result_rank)
            gene_result_p = list(gene_result_table['P-VALUE'])[0]
            data['Gene_Level_p_' +str(RUN_TYPE)].append(gene_result_p)
            gene_result_combined_score = list(gene_result_table['EXOMISER_GENE_COMBINED_SCORE'])[0]
            data['Gene_Level_combined_' +str(RUN_TYPE)].append(gene_result_combined_score)
            gene_result_variant_score = list(gene_result_table['EXOMISER_GENE_VARIANT_SCORE'])[0]
            data['Gene_Level_variant_'+str(RUN_TYPE)].append(gene_result_variant_score)
            gene_result_pheno_score = list(gene_result_table['EXOMISER_GENE_PHENO_SCORE'])[0]
            data['Gene_Level_pheno_'+str(RUN_TYPE)].append(gene_result_pheno_score)
            
            
        else:   
            print(gene, 'not prioritized for', ID)
        
            data['Gene_Level_' +str(RUN_TYPE)].append('Fail')
            data['Gene_Level_rank_' +str(RUN_TYPE)].append('N/A')
            data['Gene_Level_p_' +str(RUN_TYPE)].append('N/A')
            data['Gene_Level_combined_' +str(RUN_TYPE)].append('N/A')
            data['Gene_Level_variant_'+str(RUN_TYPE)].append('N/A')
            data['Gene_Level_pheno_'+str(RUN_TYPE)].append('N/A')

                
variant_table = pd.DataFrame(data)
print(failed)

##save results table (figure input)
variant_table.to_csv('/path/to/directory/exomiser_results_table.tsv', sep='\t', index=None)

