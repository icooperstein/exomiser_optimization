# exomiser_optimization
Example run files and scripts used in the analyses of "An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser"


## Running Exomiser / Genomiser
Installation instructions and instructions to run Exomiser and Genomiser can be found in their documentation: https://exomiser.readthedocs.io/en/latest/running.html

It is not necessary to replicate our set-up for running Exomiser or Genomiser. Example YAML files using our recommended optimized parameters can be found in the "run_exomiser/yml_files" directory for Exomiser and "run_genomiser/yml_files" for Genomiser

We provide our execution scripts to run Exomiser/Genomiser with (A) a single family or (B) cohort-level (concurrently create necessary scripts to run Exomiser/Genomiser on multiple families)


## Filtering VCFs
```
bcftools +fill-tags -O z sample.vcf.gz -- -t FORMAT/VAF,HWE | bcftools view -O z -o sample.filtered.vcf.gz -e 'FORMAT/GQ[@sample.txt] < 20 || ( GT[@sample.txt]="het" && ( FORMAT/VAF[@sample.txt] < 0.15 || FORMAT/VAF[@sample.txt] > 0.85  ) )' ```

```
* filters 0.15 ≤ VAF ≤ 0.85; GQ ≥ 20
* requires a "sample.txt" file which simply has the sample_id name (as found in VCF header)


## Subset to proband-only VCF

```
bcftools view -O u -o sample.singleton.vcf.gz -s "SAMPLE_ID" sample.filtered.vcf.gz

```

