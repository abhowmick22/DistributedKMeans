#!/usr/bin/env bash

# Output file
o=input/dnaStrandsClustPoints.csv

# Size of DNA strand
s=10

#Number of Points
n=0

#Dummy Number of Cluster
c=3

#Dummy Number of Points per cluster
p=5

		echo ********GENERATING TOTAL $n POINTS  
		python ./randomclustergen/generateDNAstrands.py -o $o -s $s -n $n -c $c -p $p

