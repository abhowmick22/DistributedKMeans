#!/usr/bin/env bash

# Output file
o=input/dnaStrandsTotPoints.csv

# Size of DNA strand
s=10

#Number of Points
n=10

#Dummy Number of Cluster
c=2

#Dummy Number of Points per cluster
p=2

		echo ********GENERATING TOTAL $n POINTS  
		python ./randomclustergen/generateDNAstrands.py -o $o -s $s -n $n -c $c -p $p

