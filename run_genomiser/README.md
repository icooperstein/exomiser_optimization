Instructions to run Exomiser and Genomiser can be found in their documentation: https://exomiser.readthedocs.io/en/latest/running.html

It is not necessary to replicate our set-up for running Exomiser/Genomiser. Example YAML files using our optimized parameters can be found in the "yml_files" directory. 

Necessary files to run Genomiser as a slurm job:
1. application.properties
2. run_exomiser.sh
3. submit_genomiser.sh 

Linux bash command: 
``` sbatch submit_genomiser.sh ID run_type ```

Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files
