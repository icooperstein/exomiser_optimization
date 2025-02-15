## Filtering VCFs
* filter_vcfs.sh
    * 0.15 ≤ VAF ≤ 0.85
    * GQ ≥ 20
* sample.txt
  * filter_vcfs.sh references a "sample.txt" file which simply has the sample_id name (as found in VCF header)

### Filtering command



```
bcftools +fill-tags -O z sample.vcf.gz -- -t FORMAT/VAF,HWE | bcftools view -O z -o sample.filtered.vcf.gz -e 'FORMAT/GQ[@sample.txt] < 20 || ( GT[@sample.txt]="het" && ( FORMAT/VAF[@sample.txt] < 0.15 || FORMAT/VAF[@sample.txt] > 0.85  ) )' ```

```
