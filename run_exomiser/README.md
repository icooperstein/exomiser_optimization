Instructions to run Exomiser and Genomiser can be found in their documentation: https://exomiser.readthedocs.io/en/latest/running.html

It is not necessary to replicate our set-up for running Exomiser/Genomiser. Example YAML files using our optimized parameters can be found in the "yml_files" directory.

We provide execution scripts to run Exomiser with a single family. Python scripts to create many scripts for cohort of families can be found in the analyses folder.

## Execution scripts for a single family
Necessary files to run Exomiser as a slurm job:
1. Proband-only or multisample VCF - provided by user
2. Pedigree (for multisample VCF) - provided by user
3. YAML file
    - [Example YAML files for Exomiser](https://github.com/icooperstein/exomiser_optimization/blob/main/run_exomiser/yml_files) 
        - optimized_exomiser.yml
            - Our final set of optimized parameters as described in publication
        - default_exomiser.yml
            - Default parameters as installed with Exomiser v14.0.0
        - no_phenotypes_exomiser.yml
            - Run Exomiser without any HPO phenotypes
        - possible_parameters_exomiser.yml
            - Commented lines describe all possible options for key parameters for use exploration
4. application.properties
5. run_exomiser.sh
6. submit_exomiser.sh

Linux bash command: ```sbatch submit_exomiser.sh ID run_type```

Replace "ID" with your sample ID and "run_type" with naming convention you have named your YML files
