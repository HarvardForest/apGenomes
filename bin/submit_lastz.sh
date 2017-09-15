#!/bin/bash
#SBATCH -n 1                  # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-02:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem-per-cpu=10000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o sub.lastz.out      # File to which STDERR will be written
#SBATCH -e sub.lastz.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=matthewklau@fas.harvard.edu # Email to which notifications will be sent

for i in $(ls ~/storage/apg/ | grep SM); 
do 
    for j in $(ls ~/storage/apg/ | grep SM); 
    do 
	sbatch lastz_batch.sh $i $j;
	sleep 5s
    done
done
