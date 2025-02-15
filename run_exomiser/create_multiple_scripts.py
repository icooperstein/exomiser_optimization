import pandas as pd
import os 

path = '/path/to/working/directory/'
phenotype_data = pd.read_csv('')

vcfs_path = "/scratch/ucgd/lustre/UCGD_Staging/UDN_Data/phase3/no_phenotype/"
peds_path = '/scratch/ucgd/lustre-labs/marth/scratch/u6013141/UDN_2024/Exomiser/updated_peds/'

RUN_TYPE = ''
id_list = [] # list of sample ids to make scripts for

failed =[]
for sample_id in id_list:
    phenotype_list1 = phenotype_data[phenotype_data['ID']==sample_id]['Terms'].item() ##replace column names to match your file with phenotype lists; comma-separated list of temrs
    #create a directory for the exomiser files for this UDN ID
    exomiser_patient_path = path + str(sample_id) + '/'
    ##also create a directory for results 
    exomiser_patient_results_path = path +  str(sample_id) + '/' + 'results/'
    if not os.path.exists(exomiser_patient_path): os.mkdir(exomiser_patient_path)
    if not os.path.exists(exomiser_patient_results_path): os.mkdir(exomiser_patient_results_path)

    ##identify vcf and ped file paths
   
    
    vcf_path = data_path  + 'vcf_name.vcf.gz'
    ped_path =  ped_dir + str(sample_id)+'.ped'

    if not os.path.exists(vcf_path): 
        if not os.path.exists(ped_path):
                failed.append(sample_id) 
        else:
            print(str(sample_id),'vcf doesnt exist but ped does')
            failed.append(sample_id)

    else:
        ##Create yml file
        yaml = open(str(exomiser_patient_path) + str(sample_id)  +'_'+str(RUN_TYPE)+'.yml', 'w')
        print('analysis:', file = yaml)
        print('  genomeAssembly: hg38', file = yaml)
        print('  vcf: ', vcf_path, sep = '', file = yaml)
        print('  ped: ', ped_path, sep = '', file = yaml)
        print('  proband: ', sample_id, sep = '', file = yaml)
        print("  hpoIds: ['" + str(phenotype_list) +"']", file = yaml)
        print('', file = yaml)
        print('  inheritanceModes: {', file = yaml)
        print('    AUTOSOMAL_DOMINANT: 0.1,', file = yaml)
        print('    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,', file = yaml)
        print('    AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,', file = yaml)
        print('    X_DOMINANT: 0.1,', file = yaml)
        print('    X_RECESSIVE_HOM_ALT: 0.1,', file = yaml)
        print('    X_RECESSIVE_COMP_HET: 2.0,', file = yaml)
        print('    MITOCHONDRIAL: 0.2', file = yaml)
        print('  }', file = yaml)
        print('', file = yaml)
        print('  analysisMode: FULL', file = yaml)
        print('', file = yaml)
        print('  frequencySources: [', file = yaml)
        print('    THOUSAND_GENOMES,', file = yaml)
        print('    TOPMED,', file = yaml)
        print('    UK10K,', file = yaml)
        print('    ESP_AFRICAN_AMERICAN, ESP_EUROPEAN_AMERICAN, ESP_ALL,', file = yaml)
        print('    EXAC_AFRICAN_INC_AFRICAN_AMERICAN, EXAC_AMERICAN,', file = yaml)
        print('    EXAC_SOUTH_ASIAN, EXAC_EAST_ASIAN,', file = yaml)
        print('    EXAC_NON_FINNISH_EUROPEAN,', file = yaml)
        print('    GNOMAD_E_AFR,', file = yaml)
        print('    GNOMAD_E_AMR,', file = yaml)
        print('    GNOMAD_E_EAS,', file = yaml)
        print('    GNOMAD_E_NFE,', file = yaml)
        print('    GNOMAD_E_SAS,', file = yaml)
        print('    GNOMAD_G_AFR,', file = yaml)
        print('    GNOMAD_G_AMR,', file = yaml)
        print('    GNOMAD_G_EAS,', file = yaml)
        print('    GNOMAD_G_NFE,', file = yaml)
        print('    GNOMAD_G_SAS', file = yaml)
        print('  ]', file = yaml)
        print('', file = yaml)
        print('  pathogenicitySources: [', file = yaml)
        # print('    CADD,', file = yaml)
        print('    REVEL,', file = yaml)
        print('    MVP,' ,file = yaml)
        print('    ALPHA_MISSENSE,', file = yaml)
        print('    SPLICE_AI', file = yaml)
        # print('    POLYPHEN,', file = yaml)
        # print('    MUTATION_TASTER,', file = yaml)
        # print('    SIFT', file = yaml)

        print('  ]', file = yaml)
        print('', file = yaml)
        
        print('  steps: [', file = yaml)
        print('    failedVariantFilter: { },', file = yaml)
        print('    variantEffectFilter: {', file = yaml)
        print('      remove: [', file = yaml)
        print('        FIVE_PRIME_UTR_EXON_VARIANT,', file = yaml)
        print('        FIVE_PRIME_UTR_INTRON_VARIANT,', file = yaml)
        print('        THREE_PRIME_UTR_EXON_VARIANT,', file = yaml)
        print('        THREE_PRIME_UTR_INTRON_VARIANT,', file = yaml)
        print('        NON_CODING_TRANSCRIPT_EXON_VARIANT,', file = yaml)
        print('        NON_CODING_TRANSCRIPT_INTRON_VARIANT,', file = yaml)
        print('        CODING_TRANSCRIPT_INTRON_VARIANT,', file = yaml)
        print('        UPSTREAM_GENE_VARIANT,', file = yaml)
        print('        DOWNSTREAM_GENE_VARIANT,', file = yaml)
        print('        INTERGENIC_VARIANT,', file = yaml)
        print('        REGULATORY_REGION_VARIANT', file = yaml)
        print('      ]', file = yaml)
        print('    },', file = yaml)
        if 'freq1' in RUN_TYPE:
            print('    frequencyFilter: {maxFrequency: 1.0},', file = yaml)
        else:
            print('    frequencyFilter: {maxFrequency: 2.0},', file = yaml)
        print('    pathogenicityFilter: {keepNonPathogenic: true},', file = yaml)
        print('    inheritanceFilter: {},', file = yaml)
        print('    omimPrioritiser: {},', file = yaml)
        if 'human' in RUN_TYPE:
            print('    hiPhivePrioritiser: {runParams: "human"},', file = yaml)
        else:
            print('    hiPhivePrioritiser: {},', file = yaml)
       # print('    hiPhivePrioritiser: {runParams: "mouse"},', file = yaml)
        #print('    phenixPrioritiser: {},', file = yaml)
        
        print('  ]', file = yaml)
        print('', file = yaml)
        print('outputOptions:', file = yaml)
        print('  outputContributingVariantsOnly: false', file = yaml)
        print('  numGenes: 0', file = yaml)
        print('  outputDirectory:  str(path), '/', str(sample_id),'/results/', sep = '', file = yaml)
        print('  outputFileName:  ', sample_id+str('_')+str(RUN_TYPE), sep = '', file = yaml)
        print('  outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]', file = yaml)
        yaml.close()

        ##Generate an exomiser script
        exFile = exomiser_patient_path +str(sample_id) + "." + str(RUN_TYPE) + '.exomiser.sh'
        exScript = open(exFile, 'w')
        print('set -eou pipefail', file = exScript)
        print('TOOLSPATH=/path/to/tool/', file = exScript) #replace with path to exomiser download directory
        print('WORKINGPATH=', str(exomiser_patient_path), sep = '', file = exScript)
        print('PROPERTIESPATH=', str(my_path), sep = '', file = exScript)
        print('STDOUT=$WORKINGPATH/', str(sample_id)+'.' +str(RUN_TYPE) + '.stdout', sep = '', file = exScript)
        print('STDERR=$WORKINGPATH/', str(sample_id) +'.'+str(RUN_TYPE) + '.stderr', sep = '', file = exScript)
        print(file = exScript)
        print('module load openjdk/17.0.1', file = exScript)
        print('java -jar $TOOLSPATH/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar \\', file = exScript)
        print('  --analysis $WORKINGPATH/', str(sample_id)+'_'+str(RUN_TYPE), '.yml \\', sep = '', file = exScript)
        print('  --spring.config.location=$PROPERTIESPATH/application.properties \\', file = exScript)
        print('  > $STDOUT \\', file = exScript)
        print('  2> $STDERR', file = exScript)
        exScript.close()

        ## Make the annotation script executable 
        makeExecutable = os.popen('chmod +x ' + str(exFile)).read()


##create execution script file with line for each exomiser run command
subFile = path + str(RUN_TYPE) + '_submission_commands.sh'
subScript = open(subFile, 'w')
for sample_id in id_list:
      print('cd', str(sample_id), ';', 
          'sbatch' ,str(path) + 'submit_exomiser.sh' ,
          str(sample_id), str(RUN_TYPE),
          '; cd ../', 
          file=subScript)


subScript.close()

