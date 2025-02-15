#! /bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=account_name ##replace with your account name
#SBATCH --partition=partition_name ##replace with your partition name

set -eou pipefail

ID=$1
RUN_TYPE=$2

bash $ID.$RUN_TYPE.genomiser.sh
