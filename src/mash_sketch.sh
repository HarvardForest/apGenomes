#!/bin/bash

### Takes input of location of main genome directory
### i.e., apg.

cd $1

for i in $(ls | grep SM); do
    mash sketch $i/scaffolds.fasta 
done

say "Hey, Mash just finished!"
