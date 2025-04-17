import pandas as pd
import altair as alt
'''
BENCHMARKING STRATEGY: Variant-level success with correct MOI (#3)
This is the most stringent success criteria.
The diagnostic variant is prioritized with the correct position, nucleotide change, and MOI. 

ie: 
chr1:123456 A>T
MOI: AR
exomiser form: 1-123456-A-T_AR

exomiser result: 1-123456-A-T-AD

This is a FAIL by this metric because the MOI is incorrect

Necessary columns in input table:
1. ID (form patientID_gene)
2. genomic_coordinates_hg38 (chr1:123456 A>T)
3. MOI (AR, AD, XR, or XD)
Note: if performing singleton analysis, all Exomiser results will be in the form 1-123456-A-T_ANY so all variants will fail unless you set MOI column to ANY 

'''

variant_table = pd.read_csv('input.txt', sep='\t')
print(len(set(variant_table['ID'])),'genes')
print(len(variant_table), 'variants')


'''Add column to variant table of variant coordinates in Exomiser's form of variant reporting'''
exomiser_variants =[]
for i, row in variant_table.iterrows():
    variant = row['genomic_coordinates_hg38']
    moi = row['MOI']
    try:
        variant_exo = variant.replace('chr', '',).replace(':','-').replace(' ','-').replace('>','-') + '_' + str(moi)
    except:
        variant_exo = 'N/A'
    exomiser_variants.append(variant_exo)
variant_table['variant_exomiser_form'] = exomiser_variants



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
        variant = row['variant_exomiser_form']

        try:
            results_file = pd.read_csv('/path/to/Exomiser/directory/exomiser_results/' + str(ID) +'/results/' + str(ID) + '_' + str(RUN_TYPE)+'.variants.tsv', sep='\t')
        except:
            print('no results file found for', str(ID))
            failed.append(ID)

        results_file = results_file[results_file['CONTRIBUTING_VARIANT']==1]
        exo_variants = list(results_file[results_file['GENE_SYMBOL']==gene]['ID'])
        exo_results = results_file[results_file['GENE_SYMBOL']==gene]

        parse_str = lambda x: x.split('_')[0]

        if gene in list(results_file['GENE_SYMBOL']):
            n+=1
            
            if variant in list(exo_results['ID']):
                data['Variant_Level_withMOI_' + str(RUN_TYPE)].append('Variant_Present_withMOI')
                variant_result = exo_results[exo_results['ID']==variant]
                rank_withMOI = min(variant_result['#RANK'])
                variant_of_interest_withMOI = variant_result[variant_result['#RANK']==rank_withMOI].head(1)
                data['Variant_Level_withMOI_rank_' + str(RUN_TYPE)].append(rank_withMOI)
                variant_p_withMOI = variant_of_interest_withMOI['P-VALUE'].item()
                data['Variant_Level_withMOI_p_' + str(RUN_TYPE)].append(variant_p_withMOI)
                variant_score_withMOI = variant_of_interest_withMOI['EXOMISER_VARIANT_SCORE'].item()
                data['Variant_Level_withMOI_variant_'+str(RUN_TYPE)].append(variant_score_withMOI)
                pheno_score_withMOI = variant_of_interest_withMOI['EXOMISER_GENE_PHENO_SCORE'].item()
                data['Variant_Level_withMOI_pheno_'+str(RUN_TYPE)].append(pheno_score_withMOI)
                combined_score_withMOI = variant_of_interest_withMOI['EXOMISER_GENE_COMBINED_SCORE'].item()
                data['Variant_Level_withMOI_combined_'+str(RUN_TYPE)].append(combined_score_withMOI)
                max_path_withMOI = variant_of_interest_withMOI['MAX_PATH'].item()
                data['Variant_Level_withMOI_maxPath_' + str(RUN_TYPE)].append(max_path_withMOI)
                max_path_source_withMOI = variant_of_interest_withMOI['MAX_PATH_SOURCE'].item()
                variant_type = variant_of_interest_withMOI['FUNCTIONAL_CLASS'].item()
                data['Variant_Type_' +str(RUN_TYPE)].append(variant_type)

                try:
                    if math.isnan(float(max_path_source_withMOI)):
                        data['Variant_Level_withMOI_maxPathSource_' + str(RUN_TYPE)].append('NONE')

                except:
                    data['Variant_Level_withMOI_maxPathSource_' + str(RUN_TYPE)].append(max_path_source_withMOI)

            else:
                print(variant, list(exo_results['ID']))
                data['Variant_Level_withMOI_' + str(RUN_TYPE)].append('Variant_Not_Present_withMOI')
                data['Variant_Level_withMOI_rank_' + str(RUN_TYPE)].append('N/A')
                data['Variant_Level_withMOI_p_' + str(RUN_TYPE)].append('N/A')
                data['Variant_Level_withMOI_variant_'+str(RUN_TYPE)].append('N/A')
                data['Variant_Level_withMOI_pheno_'+str(RUN_TYPE)].append('N/A')
                data['Variant_Level_withMOI_combined_'+str(RUN_TYPE)].append('N/A')
                data['Variant_Level_withMOI_maxPath_' + str(RUN_TYPE)].append('N/A')
                data['Variant_Level_withMOI_maxPathSource_' + str(RUN_TYPE)].append('Not_Prioritized')
                data['Variant_Type_' +str(RUN_TYPE)].append('N/A')
        else:   
        

            data['Variant_Level_withMOI_' + str(RUN_TYPE)].append('Gene_Not_in_Output')
            data['Variant_Level_withMOI_rank_' + str(RUN_TYPE)].append('N/A')
            data['Variant_Level_withMOI_p_' + str(RUN_TYPE)].append('N/A')
            data['Variant_Level_withMOI_variant_'+str(RUN_TYPE)].append('N/A')
            data['Variant_Level_withMOI_pheno_'+str(RUN_TYPE)].append('N/A')
            data['Variant_Level_withMOI_combined_'+str(RUN_TYPE)].append('N/A')
            data['Variant_Level_withMOI_maxPath_' + str(RUN_TYPE)].append('N/A')
            data['Variant_Level_withMOI_maxPathSource_' + str(RUN_TYPE)].append('Not_Prioritized')
            data['Variant_Type_' +str(RUN_TYPE)].append('N/A')
                
variant_table = pd.DataFrame(data)
print(failed)

##save results table (figure input)
variant_table.to_csv('/path/to/directory/exomiser_results_table.tsv', sep='\t', index=None)

