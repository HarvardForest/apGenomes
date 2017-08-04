#!/bin/bash

for i in $(ls | grep SM); do
    ./Mash/mash sketch $i/filtered.scaffolds.fasta 
done

say "Hey, Mash just finished!"
