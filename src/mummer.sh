#!/bin/bash
#SBATCH -n 8                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-03:00              # Runtime in D-HH:MM
#SBATCH -p general       # Partition to submit to
#SBATCH --mem-per-cpu=10000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o mummer_%j.out      # File to which STDOUT will be written
#SBATCH -e mummer_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=matthewklau@fas.harvard.edu # Email to which notifications will be sent

cd /n/home10/mklau/storage/ap_genomes/
module load MUMmer/3.23-fasrc04
mummer -mum -b -c SM-AJDMW/filtered.scaffolds.fasta SM-AZXXM/filtered.scaffolds.fasta > test.mums
