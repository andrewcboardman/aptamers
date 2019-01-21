#!/usr/bin/bash

FILE=$1
# Count kmers
for k in {3..30}; do 
	jellyfish count $FILE -m $k -s 1G -o mers_$FILE_$k
	jellyfish stats mers_${k}_0 > mers_$FILE_$k.stats
	jellyfish histo mers_${k}_0 > mers_$FILE_$k.histo
	jellyfish dump mers_${k}_0 > mers_$FILE_$k.fasta
done
