# exomiser_optimization
Scripts used in the analyses of "An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser"

* [Create cohort-level YAML and execution files](https://github.com/icooperstein/exomiser_optimization/blob/main/manuscript/analyses/create_multiple_exomiser_run_scripts.py.py)


## Filtering VCFs
Before running Exomiser or Genomiser, we recommend applying the following filters:
* GQ ≥ 20
* 0.15 ≤ VAF ≤ 0.85 for heterozygous variants
* ALT != *
```
bcftools +fill-tags -O z sample.vcf.gz -- -t FORMAT/VAF,HWE | bcftools view -O z -o sample.filtered.vcf.gz -e 'FORMAT/GQ[@sample.txt] < 20 || ( GT[@UDN789373.txt]="het" && ( FORMAT/VAF[@sample.txt] < 0.15 || FORMAT/VAF[@sample.txt] > 0.85  ) ) || ALT ="*"'

```

* requires a "sample.txt" file which simply has the sample_id name (as found in VCF header)


## Subset to proband-only VCF

```
bcftools view -O u -o sample.singleton.vcf.gz -s "SAMPLE_ID" sample.filtered.vcf.gz

```


