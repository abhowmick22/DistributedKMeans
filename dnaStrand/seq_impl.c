#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allfunc.h"

int calc_string_dist(char* point1, char* point2, int dim) {
	int retDist = 0;
	for(int i=0; i<dim; i++) {
		if(point1[i] != point2[i]){
			retDist = retDist + 1;
		}
	}
	return retDist;
}

int closest_cluster_calculator(char* point, char** cluster_centers,
								int numClusters, int dim) {
	int min_dist = calc_string_dist(point, cluster_centers[0], dim);
	int closest_center = 0;
	for(int i=1; i<numClusters; i++) {
		int temp_dist = calc_string_dist(point, cluster_centers[i], dim);
		if(temp_dist < min_dist) {
			min_dist = temp_dist;
			closest_center = i;
		}
	}
	return closest_center;
}

char** get_cluster_centers(char** points, int numPoints,
							int numClusters, int dim,
							int maxIter, float threshold) {
				
	//create initial centers = first numClusters number of
	//points
	char** cluster_centers = (char **)malloc(numClusters*sizeof(char*));
	int i, j;
	for(i=0;i<numClusters;i++) {
		cluster_centers[i] = (char *)malloc(dim*sizeof(char));
		for(j=0;j<dim;j++) {
			cluster_centers[i][j] = points[i][j];			
		}
	}
	
	//define array for each point and the cluster it belongs to
	int* pointBelongsTo = (int*)malloc(numPoints*sizeof(int));
	for(i=0;i<numPoints;i++) {
		//init with -1
		pointBelongsTo[i] = -1;
	}
	
	//allocate memory for keeping track of majority bases that belong to each position cluster
	// Use pointsInCluster and numPointsInCluster to implement the majority voting algorithm to count the majority 
	char** pointsInCluster = (char **)malloc(numClusters*sizeof(char*));
	//initialize #points count for each cluster-dim
	int** counter = (int**)malloc(numClusters*sizeof(int*));
	//initialize #points count for each cluster
	int* numPointsInCluster = (int*)malloc(numClusters*sizeof(int));
	for(i=0;i<numClusters;i++) {
		pointsInCluster[i] = (char *)malloc(dim*sizeof(char));	//because we sum all the points for each cluster
		counter[i] = (int *)malloc(dim*sizeof(int));	// counter for all positions in cluster string
	}
	
	//ratio = (# points changing cluster)/(total # points)
	float ratio = threshold+1;
	int iter = 0;	
	while(ratio > threshold && iter < maxIter) {
	
		for(i=0;i<numClusters;i++) {
			for(j=0;j<dim;j++) {
				pointsInCluster[i][j] = 'a';			// default base in string position
				counter[i][j] = 0;
			}
			numPointsInCluster[i] = 0;
		}

		//for all points, calculate the closest cluster
		int localPointsChanged = 0;
		for(i=0;i<numPoints;i++) {
			int closest_cluster = closest_cluster_calculator(points[i], cluster_centers, numClusters, dim);
			if(closest_cluster != pointBelongsTo[i]) {
				localPointsChanged++;
				pointBelongsTo[i] = closest_cluster;
			}
			
			for(j=0;j<dim;j++) {
				if(counter[closest_cluster][j] == 0){
					pointsInCluster[closest_cluster][j] = points[i][j];
					counter[closest_cluster][j]++;
				}
				else{
					if(points[i][j] == pointsInCluster[closest_cluster][j])
						counter[closest_cluster][j]++;
					else
						counter[closest_cluster][j]--;
				}
			}
			numPointsInCluster[closest_cluster]++;
						
		}
		
		//assign new cluster centers
		for(i=0; i<numClusters; i++) {
			if(numPointsInCluster[i]==0) {
				//shouldn't happen
				continue;
			}
			for(j=0; j<dim; j++) {
				cluster_centers[i][j] = pointsInCluster[i][j];
			}
		}
		
		iter++;
		ratio = ((float)totalPointsChanged)/numPoints;
	}
	
	printf("Iter=%d\n", iter);
	
	return cluster_centers;
}

int main(int argc, char* argv[]) {
	
	int numClusters = 3;
	int* numPoints = (int*)malloc(sizeof(int));
	char** init_centers = (char**)malloc(numClusters*sizeof(char*));
	int dim = 6;
	int totalNumPoints = 18;
	
	
	///afs/andrew.cmu.edu/usr11/ndhruva/public/test.txt
	char** points = readFromFileForMPI("/home/abhishek/15-640/project4/dnaStrand/test.txt", numPoints 
					   					dim, totalNumPoints, numClusters,
					   					rank, numProcs, init_centers);
	
	if(points==NULL) {
		MPI_Finalize();
		return 0;
	}
	char** cluster_centers = get_cluster_centers(points, *numPoints,
												  init_centers, numClusters,
												  dim, 1000000, 0.0);
	int i, j;
	if(rank==0) {
		printf("\n");
		for(i=0; i<numClusters; i++) {
			for(j=0; j<dim; j++) {
				printf("%c, ", init_centers[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		//print cluster centers
		for(i=0; i<numClusters; i++) {
			for(j=0; j<dim; j++) {
				printf("%c, ", cluster_centers[i][j]);
			}
			printf("\n");
		}
	}
	free(points);
	free(cluster_centers);
	free(init_centers);
	free(numPoints);
	return 0;
}
