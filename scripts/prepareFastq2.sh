#!/bin/bash
# rename fastq files and create samples.txt and possibly samples.split.txt

mkdir -p fastq
find rawdata -name '*.gz' | xargs ln -t fastq

for x in fastq/*.fastq.gz; do
    if [ `grep -o "_" <<< "$x" | wc -l` -gt 1 ]; then
        prename 's/-//g; s/(.+)_(S\d+)_L([^_])+_R([12])_([0-9]+)/$1-$2$3$5.$4/; s/(.+)_[^_]+_L([^_])+_R([12])_([0-9]+)/$1_$2$4.$3/; s/_//g; s/-/_/g' $x;
    fi;
done
paste <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/[._].*//') <('ls' fastq/*.fastq.gz | sed 's:.*/::; s/\..*//') | awk 'BEGIN { OFS = "\t" } $1 != $2 { print }' | uniq > samples.split.txt
'ls' fastq/*.fastq.gz | sed 's:.*/::; s/_.*//' | sort | uniq
