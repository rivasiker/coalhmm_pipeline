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

echo "alphabet=DNA" >> ../tmp/params.file
echo "input.sequence.multiparts=yes" >> ../tmp/params.file

echo "input.sequence.multiparts.reset=yes" >> ../tmp/params.file

echo "input.sequence.format=Fasta" >> ../tmp/params.file

echo "input.sequence.sites_to_use=all" >> ../tmp/params.file
echo "input.sequence.max_gap_allowed=50%" >> ../tmp/params.file
echo "model=GTR(a=1.0, b=1.0, c=1.0, d=1.0, e=1.0, theta=0.5, theta1 = 0.5, theta2 = 0.5)" >> ../tmp/params.file
echo "rate_distribution=Gamma(n=4, alpha=1.0)" >> ../tmp/params.file

echo "analysis=estimate" >> ../tmp/params.file

echo "optimize=yes" >> ../tmp/params.file
echo "optimization.method=fullD" >> ../tmp/params.file
echo "optimization.reparametrization=no" >> ../tmp/params.file
echo "optimization.verbose=2" >> ../tmp/params.file
echo "optimization.tolerance=0.0001" >> ../tmp/params.file
echo "optimization.max_number_f_eval=1000000" >> ../tmp/params.file
echo "optimization.max_number_iterations=2000" >> ../tmp/params.file
echo "optimization.pre=yes" >> ../tmp/params.file
echo "optimization.ignore_parameter=GTR.a,GTR.b,GTR.c,GTR.d,GTR.e,GTR.theta,GTR.theta1,GTR.theta2" >> ../tmp/params.file
echo "optimization.final=no" >> ../tmp/params.file

