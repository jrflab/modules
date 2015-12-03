#!/bin/bash
#rename fastq files and create samples.txt and possibly samples.split.txt

mkdir -p fastq
cd rawdata
find -name '*.gz' | xargs prename 's:./([A-Z]+)_([0-9]{4})/(.*)\.fastq\.gz:./$1_$2/$3$2.fastq.gz:'
find -name '*.gz' | xargs ln -t ../fastq
cd ..

for x in fastq/*.fastq.gz; do
    if [ `grep -o "_" <<< "$x" | wc -l` -gt 1 ]; then
        prename 's/-//g; s/(.+)_([NATGC]{6,8}|S\d+)_L([^_])+_R([12])_([0-9]+)/$1-$2$3$5.$4/; s/(.+)_[^_]+_L([^_])+_R([12])_([0-9]+)/$1_$2$4.$3/; s/_//g; s/-/_/g' $x;
    fi;
done
prename 's/IGO[^_]*//' fastq/*.gz
paste <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/[._].*//') <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/\..*//') | awk 'BEGIN { OFS = "\t" } $3 != $2 { print }' | uniq > samples.split.txt
'ls' fastq/*.fastq.gz | sed 's:.*/::; s/_.*//' | sort | uniq
