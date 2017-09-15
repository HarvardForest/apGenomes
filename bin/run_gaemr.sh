
#! bin/bash

## Takes input of the directory with the genomes

cd $1
for i in $(ls); do 
sbatch ~/regal/mklau/apGenomes/src/gaemr.sh $i ;
sleep 5 ;
done ;

