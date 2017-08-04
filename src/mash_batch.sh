#!/bin/bash

### use: sh mash_batch.sh <output file name> 

touch $1

for i in $(ls | grep SM); do
    for j in $(ls | grep SM); do 
	if [ "$i" != "$j" ] 
	then 
	    echo $i $j
	    say Mash starting $i $j
	    ./Mash/mash dist $i/filtered.scaffolds.fasta.msh $j//filtered.scaffolds.fasta.msh >> $1

	else
	    echo Skipping $i $j
	    say Skipping $i $j
	fi
    done
done

say "Hey, Mash is finally done!"
