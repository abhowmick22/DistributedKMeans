#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "allfunc.h"
//#include "mpi.h"

/*
 * Reads data points from file for both MPI and sequential access.
 * GP = General purpose.
 */
char** readFromFileForGP(char* filename, int dim, int numPoints) {
	printf("trying to read\n");	
	FILE* fp;
	fp = fopen(filename, "r");
	printf("%s\n", filename);
	if(fp == NULL) {
		return NULL;
	}
	
	char line[200];		//max 200 chars in line
	const char splitString[] = " ,\n";
	int numPointsSoFar = 0;
	int i, j;
	//initialize data points array
	char** points = (char**)malloc(numPoints*sizeof(char*));
	char* temp = (char*)malloc(numPoints*dim*sizeof(char));		//for contiguous 2d array
	printf("Reading points...\n");
	while(fgets(line, 200, fp) != NULL && numPointsSoFar < numPoints) {
		char* token = strtok(line, splitString);
		i=0;
		char nextDim = 'a';							// default base
		//identify starting location for this point in the 2d array
		points[numPointsSoFar] = &temp[numPointsSoFar*dim];
		while(token != NULL) {
			if(i==dim) {
				//take only the first dim dimensions per input line
				printf("Invalid input: only first %d dimensions will be considered at the %dth point.\n", dim, (numPointsSoFar+1));
				break;
			}
			if(!sscanf(token, "%c", &nextDim)) {
				//error: invalid input
				fclose(fp);
				printf("Invalid input: Encountered non-char at %dth point.\n", (numPointsSoFar+1));
				return NULL;
			}
			points[numPointsSoFar][i] = nextDim;
			token = strtok(NULL, splitString);
			i++;
		}
		if(i != dim) {
			//error: invalid input
			fclose(fp);
			printf("Invalid input: Insufficient dimensions at %dth point.\n", (numPointsSoFar+1));
			return NULL;
		}
		numPointsSoFar++;
		memset(line, 0, 200);
		//print the point read
		for(j=0;j<dim-1;j++) {
			printf("%c, ", points[numPointsSoFar-1][j]);
		}
		printf("%c\n", points[numPointsSoFar-1][dim-1]);
	}
	if(numPointsSoFar != numPoints) {
		//error: number of input points should be at least as much as mentioned by user
		printf("Invalid input: Insufficient number of points.");
		return NULL;
	}
	fclose(fp);
	printf("All points read successfully.\n");
	return points;
}

void writeToFileForGP(char* filename, char** cluster_centers, int numClusters, int dim) {
	FILE* file = fopen(filename, "w");
	int i, j;
	for(i=0; i<numClusters; i++) {
		for(j=0; j<dim-1; j++) {
			fprintf(file, "%c, ", cluster_centers[i][j]);
		}
		fprintf(file, "%c\n", cluster_centers[i][j]);
	}
	close(file);
}

void printToTerminal(char** cluster_centers, int numClusters, int dim) {
	int i, j;
	printf("-----Cluster centers-----\n");
	for(i=0; i<numClusters; i++) {
		for(j=0; j<dim-1; j++) {
			printf("%c, ", cluster_centers[i][j]);
		}
		printf("%c\n", cluster_centers[i][j]);
	}
}
