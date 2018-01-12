#!/bin/bash
#
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p general # Partition
#SBATCH --mem 10000 # Memory request
#SBATCH -t 3-00:00 # (D-HH:MM)
#SBATCH -o rep_mask.out # Standard out goes to this file
#SBATCH -e rep_mask.err # Standard err goes to this filehostname
#SBATCH --mail-type=ALL # email notification
#SBATCH --mail-user=matthewklau@fas.harvard.edu

### Run this from the directory containing a directory of un-masked genomes

module load RepeatMasker/4.0.5-fasrc04

cd $1

RepeatMasker -xsmall scaffolds.fasta
