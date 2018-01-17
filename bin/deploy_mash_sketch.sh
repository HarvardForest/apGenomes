#!/bin/bash

### use: sh mash_batch.sh <sketch directory> <output file name> 
for i in $(ls | grep SM); do
	echo Mash deploying $i
	sbatch ../bin/run_mash_sketch.sh $i
	sleep 5
done

echo "Mash sketching just finished deploying"
