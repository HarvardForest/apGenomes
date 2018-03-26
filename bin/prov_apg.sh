#!/bin/bash
#
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p general # Partition
#SBATCH --mem 10000 # Memory request
#SBATCH -t 0-00:45 # (D-HH:MM)
#SBATCH -o prov_apg.out # Standard out goes to this filehostname
#SBATCH -e prov_apg.err # Standard err goes to this filehostname

cd ../src
Rscript prov_apg.R
