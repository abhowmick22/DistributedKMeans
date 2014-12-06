#!/usr/bin/env bash
# Output file
o=input/dnaStrandsClustPoints.csv
# Size of DNA strand
s=20
#Dummy Number of Cluster
c=10
#Dummy Number of Points per cluster
p=5
python ./randomclustergen/generateDNAstrands.py -o $o -s $s -n 0 -c $c -p $p
