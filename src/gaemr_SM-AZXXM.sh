#!/bin/bash
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-01:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem=1000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o gaemr_%j.out      # File to which STDERR will be written
#SBATCH -e gaemr_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=matthewklau@fas.harvard.edu # Email to which notifications will be sent
#SBATCH -d=~/storage/ap_genomes/SM-AZXXM

module load GAEMR/1.0.1-fasrc03
GAEMR.py -c filtered.contigs.fasta --scaffolds=filtered.scaffolds.fasta --agp=filtered.agp --force
