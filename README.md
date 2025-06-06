# An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser
## Contents
- [Running Exomiser/Genomiser](#running-Exomiser-and-Genomiser)
    - [Installation](#installation)
    - [Filtering VCFs](#filtering-vcfs)
    - [Run Exomiser](#run-exomiser)
    - [Run Genomiser](#run-genomiser)
- [Analyses](#analyses)
    - [Prune HPO term lists](#prune-hpo-term-lists)
    - [Subset to proband-only VCF](#subset-to-proband-only-vcf)
    - [Cohort-level YAML and execution files](analyses/create_multiple_exomiser_run_scripts.py)
    - [Benchmarking results tables](#benchmarking-results-tables)
    - [Complete table of Exomiser/Genomiser results](#complete-table-of-results)
    - [Filtering analysis](#filtering-analysis)
    - [Remove frequently ranked genes](#remove-frequently-ranked-genes)
    - [Run time analysis](#run-time-analysis)
- [Figures](#figures)
- [Supplementary Figures](#supplementary-figures)

## Running Exomiser and Genomiser
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
    - [Example YAML files for Exomiser](run_exomiser/yml_files) 
4. [application.properties](run_exomiser/application.properties)
5. Execution scripts
    - [run_exomiser.sh](run_exomiser/run_exomiser.sh)
    - [submit_exomiser.sh](run_exomiser/submit_exomiser.sh) \
Linux bash command: ```sbatch submit_exomiser.sh ID run_type``` \
Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files
* [Create cohort-level YAML and execution files](analyses/create_multiple_exomiser_run_scripts.py)

### Run Genomiser
It is not necessary to replicate our set-up for running Genomiser. \
Necessary files to run Exomiser as a slurm job: 
1. Proband-only or multisample VCF - provided by user
2. Pedigree (for multisample VCF) - provided by user
3. YAML file
    - [Example YAML files for Genomiser](run_genomiser/yml_files) 
4. [application.properties](run_genomiser/application.properties)
5. Execution scripts
    - [run_genomiser.sh](run_genomiser/run_genomiser.sh)
    - [submit_genomiser.sh](run_genomiser/submit_genomiser.sh) \
Linux bash command: ```sbatch submit_genomiser.sh ID run_type``` \
Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files

## Analyses
### Prune HPO term lists
HPO term lists were "pruned" by removing all perinatal and prenatal HPO terms from probands comprehensive HPO term lists. A complete list of removed terms can be found in manuscript's supplementary materials or [prune_term_lists.py](analyses/prune_term_lists.py)

### Subset to proband-only VCF

```
bcftools view -O u -o sample.singleton.vcf.gz -s "SAMPLE_ID" sample.filtered.vcf.gz

```

### Benchmarking results tables
We defined three success criteria to evaluate Exomiser and Genomiser's ability to accurately prioritize diagnostic variants.
1. [Gene-level success](analyses/gene_level_benchmarking_results_table.py): The diagnostic gene is present in the prioritized Exomiser/Genomiser output, regardless of whether the variant position(s) or nucleotide change(s) contributing to the gene’s score are correct. This is the most lenient criterion.
2. [Variant-level success](variant_level_benchmarking_results_table.py): The diagnostic variant, including the correct position and nucleotide change, is prioritized, even if the mode of inheritance (MOI), such as autosomal dominant (AD) or autosomal recessive (AR), is incorrect. **This is our primary success measure and is most frequently reported in this manuscript.**
3. [Variant-level success with correct MOI](variant_level_withMOI_benchmarking_results_table.py): The diagnostic variant meets criterion 2 and has the correct MOI, making this the most stringent criterion.
For example, in compound heterozygous cases where one variant is not prioritized and the second is prioritized as an AD variant, this qualifies as a gene-level success (criterion 1), one success and one failure at the variant level (criterion 2), and two failures under variant-level success with correct MOI (criterion 3). \
\
The results of these scripts are the data inputs for all of the figures found in the manuscript, except fig 1.


### Complete table of results
This outputs a table of every variant in the results files for a cohort of patients. 
1. [Gene-Level](analyses/gene_level_complete_results_table.py): table of complete list of genes prioritized by Exomiser or Genomiser in the ".genes.tsv" results file for a cohort of patients.
2. [Variant-Level](analyses/variant_level_complete_results_table.py): table of complete list of variants prioritized by Exomiser or Genomiser in the ".variants.tsv" results file for a cohort of patients.

### Filtering analysis
[GQ_VAF_analysis.py](analyses/GQ_VAF_analysis.py) \
Analysis completed to process data for Figure 1: Evaluating VCF filtering criteria on 474 variants in combined WES and WGS cohorts.

### Remove frequently ranked genes
Supplementary Figure 14: Removal of frequently ranked genes in WGS Exomiser cohort. \
Supplementary Figure 15: Frequently ranked genes in WES Exomiser cohort. \
[remove_frequent_genes_rerank.py](analyses/remove_frequent_genes_rerank.py) \
Remove frequently prioritized genes (p≤0.3 in top 30 candidates for ≥5% of cohort) followed by subsequent reranking of diagnostic benchmarking variants. \
bash command: ```python remove_frequent_genes_rerank.py run_type WGS exomiser_results_table.tsv``` \
variables: 
1. run_type (ex: 'exomiser_default'; must match naming convention of your results files (eg run_type_sampleid.variants.tsv))
2. cohort (can be WES or WGS depending on which list of genes you want removed)
3. variants_Table (ex exomiser_results_table.tsv; created from [benchmarking results tables](#benchmarking-results-tables))

### Run time analysis
Supplementary Figure 1: Comparison of Exomiser and Genomiser performance on coding variants. \
[runtime.py](analyses/runtime.py)


## Figures
[Figure 1](figures/figure1.ipynb): Evaluating VCF filtering criteria on 474 variants in combined WES and WGS cohorts. \
[Figure 2](figures/figure2.ipynb): Evaluation of phenotype prioritization algorithms in the WGS Exomiser cohort. \
[Figure 3](figures/figure3.ipynb): Stepwise optimization process for Exomiser and Genomiser across three UDN cohorts. \
[Figure 4](figures/figure4.ipynb): Evaluation of variant pathogenicity prediction score sources in the WGS Exomiser cohort. \
[Figure 5](figures/figure5.ipynb): Frequently ranked genes in the WGS Exomiser cohort. \
[Figure 6](figures/figure6.ipynb): Parameter optimization shifts diagnostic variants into the top ten candidates in the WGS Exomiser cohort. \
Figure 7: Recommended workflow for using Exomiser and Genomiser in rare disease diagnostics. (flow chart created in Adobe Illustrator) 

## Supplementary Figures
[Supplementary Figure 1](figures/supplementary_fig1.ipynb): Comparison of Exomiser and Genomiser performance on coding variants. \
[Supplementary Figure 2](figures/supplementary_fig2.ipynb): Breakdown of three UDN cohorts used in benchmarking. \
[Supplementary Figure 3](figures/supplementary_fig3.ipynb): Summary of combined WES/WGS cohorts. \
[Supplementary Figure 4](figures/supplementary_fig4.ipynb): Genomiser performance using hiPHIVE with default versus human-only gene:phenotype associations. \
[Supplementary Figure 5](figures/supplementary_fig5.ipynb): Evaluation of individual variant pathogenicity prediction score sources in WGS Exomiser cohort. \
[Supplementary Figure 6](figures/supplementary_fig6.ipynb): WGS Exomiser diagnostic variants’ maximum pathogenicity score source broken down by variant class. \
[Supplementary Figure 7](figures/supplementary_fig7.ipynb): Evaluation of variant pathogenicity prediction score sources in WES Exomiser cohort. \
[Supplementary Figure 8](figures/supplementary_fig8.ipynb): Evaluation of variant pathogenicity prediction score sources in
WGS Genomiser cohort. \
[Supplementary Figure 9](figures/supplementary_fig9.ipynb): Exomiser performance on WGS Genomiser cohort. \
[Supplementary Figure 10](figures/supplementary_fig10.ipynb): Impact of proband phenotype quality on Exomiser performance. \
[Supplementary Figure 11](figures/supplementary_fig11.ipynb): Recovering diagnostic variants through proband-only reanalysis or manual pedigree correction. \
[Supplementary Figure 12](figures/supplementary_fig12.ipynb): Impact of family variant data and inheritance filters on Exomiser performance. \
[Supplementary Figure 13](figures/supplementary_fig13.ipynb): Impact of maximum p-value thresholds on length of candidate variant list and loss of diagnostic variants. \
[Supplementary Figure 14](figures/supplementary_fig14.ipynb): Removal of frequently ranked genes in WGS Exomiser cohort. \
[Supplementary Figure 15](figures/supplementary_fig15.ipynb): Frequently ranked genes in WES Exomiser cohort. \
[Supplementary Figure 16](figures/supplementary_fig16.ipynb): Parameter optimization shifts diagnostic variants into top ten candidates in WGS Genomiser and WES Exomiser cohorts. 