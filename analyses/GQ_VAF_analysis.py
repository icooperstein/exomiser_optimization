import pandas as pd
import gzip
from scipy.stats import rankdata
import numpy as np
import sys
run_type = sys.argv[1]

run_type= 'exomiser_default'

path = 'path/to/data/Exomiser/exomiser_results/'


def main():
    print('beginning analysis...')
    output_data = []
    output_variant_info=[]
    variant_table = pd.read_csv('input.tsv', sep='\t')
    

    patient_list = set([x.split('_')[0] for x in variant_table['ID']])

    variant_list = list(zip(variant_table['ID'],  [x.split('_')[0] for x in variant_table['variant_exomiser_form']]))


    print(len(patient_list), 'patients to analyze')
    print(len(variant_list), 'variants to analyze')
    
    
    for patient_variant in variant_list:
        sample_id = patient_variant[0].split('_')[0]
        gene = patient_variant[0].split('_')[1]
        variant = patient_variant[1]
        individual_id = get_individual_id(sample_id) 
        print(individual_id)
        '''Open and parse VCF (NO filtering)'''
        names = get_vcf_names(path +str(sample_id)+'/results/'+str(sample_id)+ '_' + str(run_type)+'.vcf.gz')
        vcf = pd.read_csv(path +str(sample_id)+ '/results/' + str(sample_id) +'_' + str(run_type) +'.vcf.gz', compression='gzip', delim_whitespace =True, comment='#', header=None, names=names, encoding='latin-1')
        
        # '''Parse VCF to get info that would be in Exomiser's variant.tsv output + GQ column'''
        header = '#RANK|ID|GENE_SYMBOL|ENTREZ_GENE_ID|MOI|P-VALUE|EXOMISER_GENE_COMBINED_SCORE|EXOMISER_GENE_PHENO_SCORE|EXOMISER_GENE_VARIANT_SCORE|EXOMISER_VARIANT_SCORE|CONTRIBUTING_VARIANT|WHITELIST_VARIANT|FUNCTIONAL_CLASS|HGVS|EXOMISER_ACMG_CLASSIFICATION|EXOMISER_ACMG_EVIDENCE|EXOMISER_ACMG_DISEASE_ID|EXOMISER_ACMG_DISEASE_NAME'.split('|')

        ##Parse INFO column for Exomiser result + GQ 
        try:
            result = parse_vcf(individual_id, vcf, header)
        except:
            'unable to parse_vcf - check sample ID'


        ###only want to rank contributing variants
        variant_tsv = pd.read_csv(path + str(sample_id) + '/results/' + str(sample_id) + '_' + str(run_type) + '.variants.tsv', sep='\t')

        checkResult(result, variant_tsv)
        result = result[result['CONTRIBUTING_VARIANT'] == '1']
        ranks = reRank(result)
        result['OG_Rank'] = ranks
        result['Variant_noMOI'] = [x.split('_')[0] for x in result['ID']]

        diag_variant_rank, diag_variant_info = find_variant_rank(variant, result, 'OG_Rank')
        #output_variant_info.append(diag_variant_info)

        output_data.append([sample_id, gene, variant, diag_variant_rank, 'OG' , 'OG', 'OG', 'OG'])
        #pd.concat(output_variant_info).to_csv(path+ str(run_type) +'_diagnostic_variant_vcf_info.tsv', sep='\t' , index=None)

        
        ###multiple restrictions
        GQs = [10, 20, 30, 40, 50, 60]
        VAFs = [[.10,.90], [.15,.85], [.20,.80], [.25,.75]]
        for min_GQ in GQs:
            print(min_GQ)
            for VAF_range in VAFs:
                print(VAF_range)
                #contributing_only
                filtered = result[result['CONTRIBUTING_VARIANT'] == '1']
                ##Remove ALT+*
                filtered = filtered[filtered['ALT']!='*']
                ##GQ
                filtered = filtered[filtered['GQ'] >= min_GQ]
                #VAF
                keep= filtered[filtered['GT']=='1/1']
                hets = filtered[filtered['GT']!='1/1']
                #print(set(hets['GT']))
                VAF_low = VAF_range[0]
                VAF_high = VAF_range[1]
                hets_filt = hets[(hets['VAF']>=VAF_low) &(hets['VAF'] <= VAF_high) ]
                hets_filt = hets_filt[hets_filt['CONTRIBUTING_VARIANT'] == '1']
                VAF_filt = pd.concat([keep, hets_filt]).sort_values(by='EXOMISER_GENE_COMBINED_SCORE',ascending=False).reset_index(drop=True)
                filtered = VAF_filt[VAF_filt['CONTRIBUTING_VARIANT'] == '1']
                #rerank
                ranks = reRank(filtered)
                filtered['Filtered_Rank'] = ranks
                diag_variant_rank,diag_variant_info = find_variant_rank(variant, filtered, 'Filtered_Rank')
                output_data.append([sample_id, gene, variant, diag_variant_rank, 'Filtered', str(min_GQ), str(VAF_range), 'no*'])

    
    
    outputDf = pd.DataFrame(output_data, columns=['ID', 'Diagnostic_Gene', 'Diagnostic_Variant','Rank', 'Class', 'GQ', 'VAF', 'Alt'])

    outfilename = 'GS_filtering_defaults_combos_' + str(run_type)+'.tsv'
    outputDf.to_csv(outfilename, sep='\t', index=None)

    print(outfilename + ' saved')



def find_gene_rank(gene, table, column):
    if gene in list(table['GENE_SYMBOL']):
        gene_rank = list(table[table['GENE_SYMBOL'] == gene][str(column)])[0]
    else:
        gene_rank = '#N/A'
    
    return gene_rank

def find_variant_rank(variant, table, column):
    if variant in list(table['Variant_noMOI']):
        variant_rank = list(table[table['Variant_noMOI'] == variant][str(column)])[0]
        variant_info = table[table['Variant_noMOI'] == variant].head(1)
    else:
        variant_rank = '#N/A'
        variant_info = pd.DataFrame()
    
    return variant_rank, variant_info


def reRank(unRankedResult):
    combined_scores = list(unRankedResult['EXOMISER_GENE_COMBINED_SCORE'])
    var_ranks = rankdata(combined_scores, method='max')
    #invert ranks for descending order
    variant_ranks  = len(combined_scores) +1 - var_ranks
    
    toRank = unRankedResult[['GENE_SYMBOL', 'MOI', 'EXOMISER_GENE_COMBINED_SCORE']]
    toRank =toRank.sort_values(by='EXOMISER_GENE_COMBINED_SCORE', ascending=False)
    toRank['GENE_MOI'] = np.array(unRankedResult['GENE_SYMBOL']) + str('_') + np.array(unRankedResult['MOI'])
    gene_list = list(toRank['GENE_MOI'])
    gene_score_data = list(zip(combined_scores, gene_list))
    running_score = 1
    running_gene = ''
    running_rank = 1
    n=0
    ranks=[]
    for score, gene in gene_score_data:
        score = float(score)
        if score == running_score:
            rank = running_rank
            if gene != running_gene:
                n +=1
                running_gene = gene                
        elif score < running_score:
            rank = running_rank + n 
            running_rank = rank 
            running_score = score
            running_gene = gene
            n=1
        else:
            print(score, running_score)
            break
        ranks.append(rank)
    # for item in list(zip(ranks, gene_list, combined_scores))[50:100]:
    #     print(item)
    print(ranks == list(unRankedResult['#RANK']))
    return ranks

    
def get_vcf_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names

def get_individual_id(sample_id):
    id_table = pd.read_csv(path+'familyIDTabletsv', sep='\t')
    individual_id = sample_id + '\n'
    return individual_id


def parse_vcf(individual_id,vcf, header):

    from collections import defaultdict
    data = defaultdict(list)
    '''Open and parse VCF (NO filtering)'''
    ##Parse INFO column for Exomiser result + GQ 
    nums = range(0, len(vcf))
    entries = []
    for i in nums:
        info = vcf.loc[i]['INFO'].split(';')
        GT = vcf.loc[i][individual_id].split(':')[0]
        GQ = int(vcf.loc[i][individual_id].split(':')[3])
        DP = int(vcf.loc[i][individual_id].split(':')[2])
        num_alt = int(vcf.loc[i][individual_id].split(':')[1].split(',')[1])
        num_ref = int(vcf.loc[i][individual_id].split(':')[1].split(',')[0])
        alt = vcf.loc[i]['ALT']
        try:
            VAF = num_alt/(num_alt+num_ref)
        except ZeroDivisionError:
            VAF=0
        for entry in info:
            if 'QD=' in entry:
                QD = entry.split('=')[1]
        for entry in info:
            if 'Exomiser' in entry:
                if '},' in entry:
                    for item in entry.split(',{'):
                        entry2 = item.replace('}','').replace('Exomiser={','').replace('{', '').split('|')
                        if len(entry2) != 18:
                            print(len(entry2),entry2)
                            print(entry)
                        j=0
                        for thing in entry2:
                            data[header[j]].append(thing)
                            j+=1
                        data['GQ'].append(GQ)
                        data['DP'].append(DP)
                        data['VAF'].append(VAF)
                        data['GT'].append(GT)
                        data['ALT'].append(alt)
                        data['QD'].append(float(QD))
                        #print(len(entry2))
        
                        
                #print(entry)
                else:
                    entry2 = entry.replace('}','').replace('Exomiser={','').split('|')
                    j=0
                    for thing in entry2:
                        data[header[j]].append(thing)
                        j+=1
                    data['GQ'].append(GQ)
                    data['DP'].append(DP)
                    data['VAF'].append(VAF)
                    data['GT'].append(GT)
                    data['ALT'].append(alt)
                    data['QD'].append(float(QD))
            
    result = pd.DataFrame(data)
    result['#RANK'] = [int(x) for x in list(result['#RANK'])]
    result = result.sort_values(by=['#RANK','MOI']).reset_index(drop=True)
    return result

def checkResult(result, variant_tsv):
    '''Make sure result matches .tsv output (only for unfiltered)....'''
    test = []
    for i, row in variant_tsv.iterrows():
        rank = row['#RANK']
        ID = row['ID']
        test.append([rank, ID])
    test2 = []
    for i, row in result.iterrows():
        rank = row['#RANK']
        ID = row['ID']
        test2.append([rank, ID])
    check = pd.DataFrame(test, columns=['#RANK', 'Variant']).sort_values(by='Variant').reset_index(drop=True) == pd.DataFrame(test2, columns=['#RANK', 'Variant']).sort_values(by='Variant').reset_index(drop=True)
    if False in list(check['#RANK']) or False in list(check['Variant']):
        return 'bad'
    else:
        return 'Check passed'






if __name__ == "__main__":
    main()