#!/bin/bash
#SBATCH -n 8                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-05:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem-per-cpu=10000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o b.lastz_test_1.out      # File to which STDERR will be written
#SBATCH -e b.lastz_test_1.err      # File to which STDERR will be written
#SBATCH --mail-type=END         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=matthewklau@fas.harvard.edu # Email to which notifications will be sent

~/lastz-distrib/bin/lastz ~/storage/apg/SM-AZXXM/SM-AZXXM.fna ~/storage/apg/SM-AZXXM/SM-AZXXM.fna --notransition --step=1 --nogapped --format=rdotplot > ~/storage/apg/lastzr/azxxm_vs_azxxm1.txt
