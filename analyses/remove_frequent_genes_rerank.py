import pandas as pd
import gzip
from scipy.stats import rankdata
import numpy as np
import sys
run_type = sys.argv[1] ## "ex: exomiser_default"
cohort = sys.argv[2] ## WES or WGS
variants_table = sys.argv[3] ##ex: exomiser_results_table.tsv


path = '/path/to/data/Exomiser/exomiser_results/'


def main():
    print('beginning analysis...')
    outputData = []
    output_variant_info=[]
    patient_list = set([x.split('_')[0] for x in variants_table['ID']])
    #variant_list = list(zip(variants_table['ID'],  [x.split('_')[0] for x in variants_table['variant_exomiser_form']]))
    variant_list = list(zip(variants_table['ID'], variants_table['genomic_coordinates_hg38']))
    print(len(patient_list), 'patients to analyze')
    print(len(variant_list), 'variants to analyze')
    
    id_table = pd.read_csv(path+'4.8.24_update/1723_UDN_inGenomeJC_familyIDTable_4.8.24.tsv', sep='\t')
    id_table = pd.read_csv('/scratch/ucgd/lustre-labs/marth/scratch/u6013141/UDN_2024/Exomiser_WES/229_variants_phenotips_data_10.01.24_HGVS.tsv', sep='\t',encoding='unicode_escape')

    for patient_variant in variant_list:
        sample_id = patient_variant[0].split('_')[0]
        gene = patient_variant[0].split('_')[1]
        variant = patient_variant[1]
        individual_id = get_individual_id(sample_id)
        print(individual_id)
        '''Open and parse VCF (NO filtering)'''
        names = get_vcf_names(path  +str(sample_id)+'/results/'+str(sample_id)+ '_' + str(run_type)+'.vcf.gz')
        vcf = pd.read_csv(path + str(sample_id)+ '/results/' + str(sample_id) +'_' + str(run_type) +'.vcf.gz', compression='gzip', delim_whitespace =True, comment='#', header=None, names=names, encoding='latin-1')
        
        # '''Parse VCF to get info that would be in Exomiser's variant.tsv output + GQ column'''
        header = '#RANK|ID|GENE_SYMBOL|ENTREZ_GENE_ID|MOI|P-VALUE|EXOMISER_GENE_COMBINED_SCORE|EXOMISER_GENE_PHENO_SCORE|EXOMISER_GENE_VARIANT_SCORE|EXOMISER_VARIANT_SCORE|CONTRIBUTING_VARIANT|WHITELIST_VARIANT|FUNCTIONAL_CLASS|HGVS|EXOMISER_ACMG_CLASSIFICATION|EXOMISER_ACMG_EVIDENCE|EXOMISER_ACMG_DISEASE_ID|EXOMISER_ACMG_DISEASE_NAME'.split('|')


        ##Parse INFO column for Exomiser result + GQ 
        try:
            result = parse_vcf(individual_id, vcf, header)
        except:
            individual_id = list(id_table[id_table['sample_id']  == sample_id]['sample_id'])[0]
            result = parse_vcf(individual_id, vcf, header)
        variantTsv = pd.read_csv(path + str(sample_id) + '/results/' + str(sample_id) + '_' + str(run_type) + '.variants.tsv', sep='\t')
        checkResult(result, variantTsv)

        result = result[result['CONTRIBUTING_VARIANT'] == '1']
        ranks = reRank(result)
        result['OG_Rank'] = ranks
        result['Variant_noMOI'] = [x.split('_')[0] for x in result['ID']]

        diag_variant_rank, diag_variant_info = find_variant_rank(variant, result, 'OG_Rank')
        #output_variant_info.append(diag_variant_info)

        outputData.append([sample_id, gene, variant, diag_variant_rank, 'OG' ])
        #pd.concat(output_variant_info).to_csv('/scratch/ucgd/lustre-labs/marth/scratch/u6013141/02_exomiser_manuscript/genomiser/' + str(run_type) +'_diagnostic_variant_vcf_info.tsv', sep='\t' , index=None)
        cohort='WES'
        ## FILTER OUT COMMON GENES###
        '''WGS Exomiser cohort: 89 genes (p ≤ 0.3; present in the top 30 candidates for 5% of cohort)'''
        '''WES Exomiser cohort: 99 genes (p ≤ 0.3; present in the top 30 candidates for 5% of cohort)'''

        if cohort == 'WGS' or cohort =='GS':
            gene_list = ['PRAMEF10', 'JMY', 'OR2T35', 'KRT18', 'MUC3A', 'OR8U1', 'SLC25A5', 'ADPRHL1', 'HLA-DRB1', 'ZNF880', 'ERICH2', 'ASAH2B', 'MAGEC1', 'PABPC3', 'KRT6B', 'USP17L18', 'TTN', 'RBMX', 'DRD4', 'ANAPC1', 'SKA3', 'GOLGA6L22', 'MUC4', 'USP17L17', 'SPTBN4', 'KCNJ18', 'FOXD4L5', 'GOLGA6L2', 'ANKRD36', 'FAM86B2', 'FAM86B1', 'TRIM64B', 'TUBB8B', 'NBPF8', 'VCX', 'MUC6', 'WASHC1', 'TYRO3', 'DUX4', 'DCAF15', 'SIRPA', 'VPS13B', 'NEB', 'DSPP', 'SHROOM4', 'KRTAP10-6', 'TPSB2', 'PLIN4', 'ZNG1C', 'IL32', 'ANKRD36C', 'FOXD4L3', 'MTCH2', 'SYNE1', 'OR8G1', 'GOLGA6L9', 'NOC4L', 'GOLGA6L10', 'USP17L10', 'CNTNAP3B', 'NXNL1', 'DGKQ', 'AKR1C8', 'HLA-DRB5', 'KRTAP5-8', 'LILRB1', 'TDG', 'OBSCN', 'HRCT1', 'PLEC', 'INF2', 'KMT2C', 'ZNF717', 'FRG2C', 'PRSS2', 'TRRAP', 'PCDHA4', 'FANCD2', 'CEL', 'KMT2B', 'USP17L11', 'GOLGA6L6', 'PRAMEF33', 'FOXD4L4', 'ARHGEF5', 'TSHZ1']

        elif cohort =='WES' or cohort=='ES':
            gene_list = ['MUC4', 'GJC2', 'TDG', 'NFIA', 'SREBF1', 'TTN', 'MYH6', 'AP3S1', 'ATAT1', 'HLA-DRB1', 'PRKRA', 'ANAPC1', 'MRPL4', 'ASXL3', 'PRPF40B', 'TYRO3', 'TRA2A', 'VPS13B', 'SKA3', 'CCDC124', 'WIPF3', 'BICRA', 'MED13L', 'ATP2B3', 'NTN1', 'CEL', 'NPIPB6', 'FOXD4L5', 'TRRAP', 'CHD5', 'RBMX', 'GOLGA6L2', 'NDUFB11', 'KMT2E', 'HLA-DRB5', 'KCNJ18', 'PLA2G4D', 'RUVBL2', 'ADM', 'OR2T35', 'TUBB8B', 'PRSS1', 'KMT2C', 'FANCD2', 'HLA-B', 'MAGEC1', 'SEC63', 'SHROOM4', 'KRTAP10-6', 'SH2B1', 'USP9X', 'RPL10', 'SIRPA', 'RIN3', 'SPTBN4', 'ZNF880', 'ZXDA', 'FAM86B1', 'ZNG1C', 'ASAH2B', 'MUC3A', 'HLA-A', 'ZNF717', 'SYNE1', 'CP', 'MUC6', 'ADAMTS7', 'PLEC', 'GOLGA6L6', 'DCAF15', 'NEB', 'PRAMEF33', 'RERE', 'KMT2B', 'IRF2BPL', 'SCN4A', 'CACNA1H', 'GLE1', 'INF2', 'RP1L1', 'OBSCN', 'VILL', 'IL32', 'VCX', 'FOXD4L3', 'HLA-DQB1', 'CACNA1I', 'FHOD3', 'SLC34A1', 'MAPK12', 'SYNGAP1', 'KRTAP5-5', 'HSPG2', 'CENPB', 'KRTAP4-8', 'TNXB', 'PCDHA4', 'HIP1', 'PRB3']
        
        print(len(gene_list), 'genes to remove based on', cohort, 'cohort')

        parse_str = lambda x: float(x)
        result['P-VALUE'] = result['P-VALUE'].map(parse_str)
        filtered=result[result['CONTRIBUTING_VARIANT']=='1']
        filtered = filtered[filtered['P-VALUE'] <= 0.3]
        filtered = filtered[~filtered['GENE_SYMBOL'].isin(gene_list)]
        ranks = reRank(filtered)
        filtered['gene_p_filt_Rank'] = ranks
        
        diag_variant_rank, diag_variant_info  = find_variant_rank(variant, filtered, 'gene_p_filt_Rank')
        outputData.append([sample_id, gene,variant, diag_variant_rank, 'gene_p_filt'])
        
        df = pd.DataFrame(outputData, columns=['ID', 'Diagnostic_Gene', 'Diagnostic_Variant','Rank', 'Class'])
        print(df[df['ID'] == sample_id])
    
    outputDf = pd.DataFrame(outputData, columns=['ID', 'Diagnostic_Gene', 'Diagnostic_Variant','Rank', 'Class'])

    outfilename = path+ '5percent_common_genes_p0.3_removed_' + str(run_type)+'.tsv'
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

def checkResult(result, variantTSV):
    '''Make sure result matches .tsv output (only for unfiltered)....'''
    test = []
    for i, row in variantTSV.iterrows():
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