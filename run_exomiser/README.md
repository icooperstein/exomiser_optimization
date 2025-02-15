Instructions to run Exomiser and Genomiser can be found in their documentation: https://exomiser.readthedocs.io/en/latest/running.html

It is not necessary to replicate our set-up for running Exomiser/Genomiser. Example YAML files using our optimized parameters can be found in the "yml_files" directory.

We provide execution scripts to run Exomiser with (A) a single family or (B) create many scripts for cohort of families

## (A) Execution scripts for a single family
Necessary files to run Exomiser as a slurm job:

1. application.properties
2. run_exomiser.sh
3. submit_exomiser.sh

Linux bash command: ```sbatch submit_exomiser.sh ID run_type```
Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files


## (B) Cohort-level script preparation and execution files
- this adaption requires '''submit_exomiser.sh'''
#### create_multiple_scripts.py
Necessary inputs: 
1. List of sample ids for which you wish to run Exomiser on
2. table of comma-separated phenotype lists for each sample id


This script will write 3 files to designated folder for each sample id:
1. YAML file: /exomiser/patient/path/sample_id_runtype.yml
2. /exomiser/patient/path/sample_id.runtype.exomiser.sh (equivalent to run_exomiser.sh)
3. /working/directory/runtype_submission_commands.sh ## one line per sample_id slurm job submission
