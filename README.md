# An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser
### Contents
- [Running Exomiser/Genomiser](https://github.com/icooperstein/exomiser_optimization#Installation-and-running-Exomiser-and-Genomiser)
    - [Installation](https://github.com/icooperstein/exomiser_optimization#installation)
    - [Filtering VCFs](https://github.com/icooperstein/exomiser_optimization#filtering-vcfs)
    - [Run Exomiser](https://github.com/icooperstein/exomiser_optimization#run-exomiser)
    - [Run Genomiser](https://github.com/icooperstein/exomiser_optimization#run-genomiser)
- [Figures](https://github.com/icooperstein/exomiser_optimization#figures)

- [Create cohort-level YAML and execution files](https://github.com/icooperstein/exomiser_optimization/blob/main/manuscript/analyses/create_multiple_exomiser_run_scripts.py.py)

## Figures
* jupyter notebooks used for figure generation in manuscript as well as PDFs of final figures


## Running Exomiser/Genomiser
### Installation
Installation instructions and instructions to run Exomiser and Genomiser can be found in their documentation: https://exomiser.readthedocs.io/en/latest/running.html
### Filtering VCFs
Before running Exomiser or Genomiser, we recommend applying the following filters to your VCF to remove potential false positive variants:
* GQ ≥ 20
* 0.15 ≤ VAF ≤ 0.85 for heterozygous variants
* ALT != *
* requirements: a "sample.txt" file which simply has the sample_id name (as found in VCF header)

```
bcftools +fill-tags -O z sample.vcf.gz -- -t FORMAT/VAF,HWE | bcftools view -O z -o sample.filtered.vcf.gz -e 'FORMAT/GQ[@sample.txt] < 20 || ( GT[@UDN789373.txt]="het" && ( FORMAT/VAF[@sample.txt] < 0.15 || FORMAT/VAF[@sample.txt] > 0.85  ) ) || ALT ="*"'

```

### Run Exomiser
It is not necessary to replicate our set-up for running Exomiser or Genomiser. \
Necessary files to run Exomiser as a slurm job: 
1. Proband-only or multisample VCF - provided by user
2. Pedigree (for multisample VCF) - provided by user
3. YAML file
    - [Example YAML files for Exomiser](https://github.com/icooperstein/exomiser_optimization/blob/main/run_exomiser/yml_files) 

4. [application.properties](https://github.com/icooperstein/exomiser_optimization/blob/main/run_exomiser/application.properties)
5. Execution scripts
    - [run_exomiser.sh](https://github.com/icooperstein/exomiser_optimization/blob/main/run_exomiser/run_exomiser.sh)
    - [submit_exomiser.sh](https://github.com/icooperstein/exomiser_optimization/blob/main/run_exomiser/submit_exomiser.sh)
Linux bash command: ```sbatch submit_exomiser.sh ID run_type```
Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files


### Run Genomiser
It is not necessary to replicate our set-up for running Genomiser. \
Necessary files to run Exomiser as a slurm job: 
1. Proband-only or multisample VCF - provided by user
2. Pedigree (for multisample VCF) - provided by user
3. YAML file
    - [Example YAML files for Genomiser](https://github.com/icooperstein/exomiser_optimization/blob/main/run_genomiser/yml_files) 

4. [application.properties](https://github.com/icooperstein/exomiser_optimization/blob/main/un_genomiser/application.properties)
5. Execution scripts
    - [run_genomiser.sh](https://github.com/icooperstein/exomiser_optimization/blob/main/run_genomiser/run_genomiser.sh)
    - [submit_genomiser.sh](https://github.com/icooperstein/exomiser_optimization/blob/main/un_genomiser/submit_genomiser.sh)
Linux bash command: ```sbatch submit_genomiser.sh ID run_type```
Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files

## Subset to proband-only VCF

```
bcftools view -O u -o sample.singleton.vcf.gz -s "SAMPLE_ID" sample.filtered.vcf.gz

```


