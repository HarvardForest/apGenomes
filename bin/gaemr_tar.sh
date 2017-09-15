#! bin/bash

### Input is the path to the directory with the genomes

cd $1 
for i in $(ls | grep SM-A); do 
cd $i ;
cp gaemr gaemr_$i ;
cd gaemr_$i ;
cd work ;
rm *.fasta ;
cd .. ;
cd .. ;
tar -zcvf gaemr_$i.tar.gz gaemr_$i ;
cd .. ;
done ;
