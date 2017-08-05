#!/bin/bash

### use: sh mash_batch.sh <output file name> 

cd $2
touch $1

for i in $(ls | grep SM); do
    for j in $(ls | grep SM); do 
	if [ "$i" != "$j" ] 
	then 
	    echo $i $j
	    say Mash starting $i $j
	    mash dist $i/scaffolds.fasta.msh $j//scaffolds.fasta.msh >> $1

	else
	    echo Skipping $i $j
	    say Skipping $i $j
	fi
    done
done

say "Hey, Mash is finally done!"
