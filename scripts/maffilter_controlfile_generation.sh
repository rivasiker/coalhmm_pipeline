#!/bin/bash

### This script create a control file for maffiler ###
### It takes five arguments: a maf file to be filtered and four species names

maffile=$1
species1=$2
species2=$3
species3=$4
species4=$5

echo "input.file=$1" > ../tmp/control_file
echo "input.file.compression=gzip" >> ../tmp/control_file
echo "input.format=Maf" >> ../tmp/control_file
echo "output.log=../tmp/maf_filtering.log" >> ../tmp/control_file
echo -e "maf.filter=Subset(species=($2, $3, $4, $5), strict=yes, keep=no, remove_duplicates=yes), Merge(species=($2, $3, $4, $5), dist_max=100, rename_chimeric_chromosomes=yes), XFullGap(species=($2, $3, $4, $5)), MinBlockLength(min_length=500), WindowSplit(preferred_size=100000, align=adjust, keep_small_blocks=yes), Output(file=../tmp/filtered_maf, compression=none, mask=no)" >> ../tmp/control_file

