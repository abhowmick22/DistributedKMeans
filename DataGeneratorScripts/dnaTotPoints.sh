#!/usr/bin/env bash

# Output file
o=input/dnaStrandsTotPoints.csv

# Size of DNA strand
s=5

#Number of Points
n=10

#Dummy Number of Cluster
c=0

#Dummy Number of Points per cluster
p=0

		echo ********GENERATING TOTAL $n POINTS  
		python ./randomclustergen/generateDNAstrands.py -o $o -s $s -n $n -c $c -p $p

