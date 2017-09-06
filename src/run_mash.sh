#!/bin/bash
#SBATCH -n 8                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-24:00              # Runtime in D-HH:MM
#SBATCH -p serial_requeue       # Partition to submit to
#SBATCH --mem-per-cpu=10000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o mash.out      # File to which STDERR will be written
#SBATCH -e mash.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=matthewklau@fas.harvard.edu # Email to which notifications will be sent

cd $1

### Run the sketch first to save time

for i in $(ls | grep SM); do
    ./Mash/mash sketch $i/filtered.scaffolds.fasta 
done

say "Hey, Mash sketiching just finished!"

### Calculate the Mash distances

for i in $(ls | grep SM); do
    for j in $(ls | grep SM); do 
	    echo $i $j
	    say Mash starting $i $j
	    ./Mash/mash dist $i/filtered.scaffolds.fasta.msh $j/filtered.scaffolds.fasta.msh >> $1
done

say "Hey, Mash is finally done!"
