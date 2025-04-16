# An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser
### Contents
- [Running Exomiser/Genomiser](https://github.com/icooperstein/exomiser_optimization#Installation-and-running-Exomiser-and-Genomiser)
    - [Installation](https://github.com/icooperstein/exomiser_optimization#installation)
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

### Run Exomiser/Genomiser
It is not necessary to replicate our set-up for running Exomiser or Genomiser. \
[Example YAML files for Exomiser](https://github.com/icooperstein/exomiser_optimization/blob/main/run_exomiser/yml_files) \
    - optimized_exomiser.yml
        - Our final set of optimized parameters as described in publication
    - default_exomiser.yml
        - Default parameters as installed with Exomiser v14.0.0
    - no_phenotypes_exomiser.yml
        - Run Exomiser without and HPO phenotypes
    - possible_parameters_exomiser.ymlhttps://github.com/icooperstein/exomiser_optimization/blob/main/run_exomiser/yml_files/possible_parameters_exomiser.yml
        - Commented lines describe all possible options for key parameters for use exploration

[Example YAML files for Genomiser](https://github.com/icooperstein/exomiser_optimization/blob/main/run_genomiser/yml_files)


## Subset to proband-only VCF

```
bcftools view -O u -o sample.singleton.vcf.gz -s "SAMPLE_ID" sample.filtered.vcf.gz

```


