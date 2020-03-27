#!/bin/bash

### This script create a file containing the parameters for coalhmm
### It takes five arguments: a maf file to be filtered and four species names

species1=$1
species2=$2
species3=$3
species4=$4

echo "coalmethod=ILS(\\" > ../tmp/params.file
echo "    implementation=09,\\" >> ../tmp/params.file
echo "    nbSpecies=3,\\" >> ../tmp/params.file
echo "    species1=$1,\\" >> ../tmp/params.file
echo "    species2=$2,\\" >> ../tmp/params.file
echo "    species3=$3,\\" >> ../tmp/params.file
echo "    outgroup=$4,\\" >> ../tmp/params.file
echo "    tau1=0.004,\\" >> ../tmp/params.file
echo "    tau2=0.0015,\\" >> ../tmp/params.file
echo "    c2=0.010,\\" >> ../tmp/params.file
echo "    theta1=0.002,\\" >> ../tmp/params.file
echo "    theta2=0.002,\\" >> ../tmp/params.file
echo "    median=no,\\" >> ../tmp/params.file
echo "    rho=0.2,\\" >> ../tmp/params.file
echo "    tau.min = 0.0001,\\" >> ../tmp/params.file
echo "    theta.min = 0.0001,\\" >> ../tmp/params.file
echo "    rho.min = 0.0001,\\" >> ../tmp/params.file
echo "    rho.max = 1000\\" >> ../tmp/params.file
echo "  )" >> ../tmp/params.file


