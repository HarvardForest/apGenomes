#!/bin/bash
#SBATCH -n 8                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-02:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem-per-cpu=10000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o b.lastz%A_%a.out      # File to which STDERR will be written
#SBATCH -e b.lastz%A_%a.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=matthewklau@fas.harvard.edu # Email to which notifications will be sent

~/lastz-distrib/bin/lastz ~/storage/apg/$1/$1.fna ~/storage/apg/$2/$2.fna --notransition --step=20 --nogapped --format=maf > ~/storage/apg/$1_vs_$2.maf


