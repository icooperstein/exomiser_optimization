Instructions to run Exomiser and Genomiser can be found in their documentation: https://exomiser.readthedocs.io/en/latest/running.html

It is not necessary to replicate our set-up for running Exomiser/Genomiser. Example YAML files using our optimized parameters can be found in the "yml_files" directory.

We provide execution scripts to run Exomiser with a single family. Python scripts to create many scripts for cohort of families can be found in the analyses folder.

## Execution scripts for a single family
Necessary files to run Exomiser as a slurm job:

1. application.properties
2. run_exomiser.sh
3. submit_exomiser.sh

Linux bash command: ```sbatch submit_exomiser.sh ID run_type```

Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files
