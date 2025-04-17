import pandas as pd
import altair as alt
'''
BENCHMARKING STRATEGY: Variant-level success (#2)
This is our primary success measure and is most frequently reported in the manuscript.
The diagnostic variant, including the correct position and nucleotide change, is prioritized, even if the MOI is incorrect.
ie: 
chr1:123456 A>T
MOI: AR
exomiser form: 1-123456-A-T_AR

exomiser result: 1-123456-A-T-AD

This is a SUCCESS by this metric; despite the MOI being incorrect 

Necessary columns in input table:
1. ID (form patientID_gene)
2. genomic_coordinates_hg38 (chr1:123456 A>T)
3. MOI (AR, AD, XR, or XD)

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
        variant_noMOI = variant.split('_')[0]

        try:
            results_file = pd.read_csv('/path/to/Exomiser/directory/exomiser_results/' + str(ID) +'/results/' + str(ID) + '_' + str(RUN_TYPE)+'.variants.tsv', sep='\t')
        except:
            print('no results file found for', str(ID))
            failed.append(ID)

        results_file = results_file[results_file['CONTRIBUTING_VARIANT']==1]
        exo_variants = list(results_file[results_file['GENE_SYMBOL']==gene]['ID'])
        exo_results = results_file[results_file['GENE_SYMBOL']==gene]

        parse_str = lambda x: x.split('_')[0]
        exo_results.loc[:,'ID2'] = exo_results.loc[:,'ID'].map(parse_str) ##for no MOI

        if gene in list(results_file['GENE_SYMBOL']):
            n+=1
            
            if variant_noMOI in list(exo_results['ID2']):
                data['Variant_Level_noMOI_' + str(RUN_TYPE)].append('Variant_Present_noMOI')
                variant_result = exo_results[exo_results['ID2']==variant_noMOI]
                rank_noMOI = min(variant_result['#RANK'])
                variant_of_interest_noMOI = variant_result[variant_result['#RANK']==rank_noMOI].head(1)
                data['Variant_Level_noMOI_rank_' + str(RUN_TYPE)].append(rank_noMOI)
                variant_p_noMOI = variant_of_interest_noMOI['P-VALUE'].item()
                data['Variant_Level_noMOI_p_' + str(RUN_TYPE)].append(variant_p_noMOI)
                variant_score_noMOI = variant_of_interest_noMOI['EXOMISER_VARIANT_SCORE'].item()
                data['Variant_Level_noMOI_variant_'+str(RUN_TYPE)].append(variant_score_noMOI)
                pheno_score_noMOI = variant_of_interest_noMOI['EXOMISER_GENE_PHENO_SCORE'].item()
                data['Variant_Level_noMOI_pheno_'+str(RUN_TYPE)].append(pheno_score_noMOI)
                combined_score_noMOI = variant_of_interest_noMOI['EXOMISER_GENE_COMBINED_SCORE'].item()
                data['Variant_Level_noMOI_combined_'+str(RUN_TYPE)].append(combined_score_noMOI)
                max_path_noMOI = variant_of_interest_noMOI['MAX_PATH'].item()
                data['Variant_Level_noMOI_maxPath_' + str(RUN_TYPE)].append(max_path_noMOI)
                max_path_source_noMOI = variant_of_interest_noMOI['MAX_PATH_SOURCE'].item()
                variant_type = variant_of_interest_noMOI['FUNCTIONAL_CLASS'].item()
                data['Variant_Type_' +str(RUN_TYPE)].append(variant_type)

                try:
                    if math.isnan(float(max_path_source_noMOI)):
                        data['Variant_Level_noMOI_maxPathSource_' + str(RUN_TYPE)].append('NONE')

                except:
                    data['Variant_Level_noMOI_maxPathSource_' + str(RUN_TYPE)].append(max_path_source_noMOI)

            else:
                print(variant_noMOI, list(exo_results['ID2']))

                data['Variant_Level_noMOI_' + str(RUN_TYPE)].append('Variant_Not_Present_noMOI')
                data['Variant_Level_noMOI_rank_' + str(RUN_TYPE)].append('N/A')
                data['Variant_Level_noMOI_p_' + str(RUN_TYPE)].append('N/A')
                data['Variant_Level_noMOI_variant_'+str(RUN_TYPE)].append('N/A')
                data['Variant_Level_noMOI_pheno_'+str(RUN_TYPE)].append('N/A')
                data['Variant_Level_noMOI_combined_'+str(RUN_TYPE)].append('N/A')
                data['Variant_Level_noMOI_maxPath_' + str(RUN_TYPE)].append('N/A')
                data['Variant_Level_noMOI_maxPathSource_' + str(RUN_TYPE)].append('Not_Prioritized')
                data['Variant_Type_' +str(RUN_TYPE)].append('N/A')
        else:   
        

            data['Variant_Level_noMOI_' + str(RUN_TYPE)].append('Gene_Not_in_Output')
            data['Variant_Level_noMOI_rank_' + str(RUN_TYPE)].append('N/A')
            data['Variant_Level_noMOI_p_' + str(RUN_TYPE)].append('N/A')
            data['Variant_Level_noMOI_variant_'+str(RUN_TYPE)].append('N/A')
            data['Variant_Level_noMOI_pheno_'+str(RUN_TYPE)].append('N/A')
            data['Variant_Level_noMOI_combined_'+str(RUN_TYPE)].append('N/A')
            data['Variant_Level_noMOI_maxPath_' + str(RUN_TYPE)].append('N/A')
            data['Variant_Level_noMOI_maxPathSource_' + str(RUN_TYPE)].append('Not_Prioritized')
            data['Variant_Type_' +str(RUN_TYPE)].append('N/A')
                
variant_table = pd.DataFrame(data)
print(failed)

##save results table (figure input)
variant_table.to_csv('/path/to/directory/exomiser_results_table.tsv', sep='\t', index=None)