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
echo "maf.filter=Subset(species=($2, $3, $4, $5), strict=yes, keep=no, remove_duplicates=yes), \\" >> ../tmp/control_file
echo "Merge(species=($2, $3, $4, $5), dist_max=100, rename_chimeric_chromosomes=yes), \\" >> ../tmp/control_file
echo "XFullGap(species=($2, $3, $4, $5)), \\" >> ../tmp/control_file
echo "MinBlockLength(min_length=500), \\" >> ../tmp/control_file
echo "AlnFilter(species=($2, $3, $4, $5), window.size=500, window.step=500, max.gap=0.2, max.ent=0.2, missing_as_gap=yes, relative=yes, file=none), \\" >> ../tmp/control_file  
echo "MinBlockLength(min_length=500), \\" >> ../tmp/control_file
echo "WindowSplit(preferred_size=100000, align=adjust, keep_small_blocks=yes), \\" >> ../tmp/control_file
echo "Output(file=../tmp/filtered.maf, compression=none, mask=no)" >> ../tmp/control_file
