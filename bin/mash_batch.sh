#!/bin/bash

### use: sh mash_batch.sh <sketch directory> <output file name> 

cd $1
touch $2

for i in $(ls | grep SM); do
    for j in $(ls | grep SM); do 
	echo $i $j
	say Mash starting $i $j
	mash dist $i/scaffolds.fasta.msh $j/scaffolds.fasta.msh >> $2
    done
done

say "Mash just finished"
