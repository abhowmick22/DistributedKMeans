#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allfunc.h"

double calc_squared_dist(double* point1, double* point2, int dim) {
	double retDist = 0.0f;
	int i=0;
	for(i=0; i<dim; i++) {
		retDist = retDist + pow((point2[i]-point1[i]), 2);
	}
	return retDist;
}

int closest_cluster_calculator(double* point, double** cluster_centers,
								int numClusters, int dim) {
	double min_dist = calc_squared_dist(point, cluster_centers[0], dim);
	int closest_center = 0;
	int i=1;
	for(; i<numClusters; i++) {
		double temp_dist = calc_squared_dist(point, cluster_centers[i], dim);
		if(temp_dist < min_dist) {
			min_dist = temp_dist;
			closest_center = i;
		}
	}
	return closest_center;
}

double** get_cluster_centers_seq(double** points, int numPoints,
							int numClusters, int dim,
							int maxIter, int totalNumPoints,
							double threshold) {
				
	//create initial centers = first numClusters number of
	//points
	double** cluster_centers = (double **)malloc(numClusters*sizeof(double*));
	int i, j;
	
	//assign first numClusters points as initial centers
//	for(i=0;i<numClusters;i++) {
//		cluster_centers[i] = (double *)malloc(dim*sizeof(double));
//		for(j=0;j<dim;j++) {
//			cluster_centers[i][j] = points[i][j];			
//		}
//	}
	
	//random initial centers
	int* prevCenters = (int*)malloc(numClusters*sizeof(int));
	for(i=0; i<numClusters; i++) {
		prevCenters[i] = -1;
	}	
	for(i=0;i<numClusters;) {
		srand(time(NULL));	//will introduce delay, but will be more randomized
		int r = rand() % totalNumPoints;
		//check if the same has been encountered before
		for(j=0;j<i;j++) {
			if(prevCenters[j]==r) {
				break;
			}
		}
		if(j!=i) {
			continue;
		}
		cluster_centers[i] = (double *)malloc(dim*sizeof(double));
		prevCenters[i] = r;
		for(j=0;j<dim;j++) {
			cluster_centers[i][j] = points[r][j];
		}
		i++;
	}
	free(prevCenters);
	
	//define array for each point and the cluster it belongs to
	int* pointBelongsTo = (int*)malloc(numPoints*sizeof(int));
	for(i=0;i<numPoints;i++) {
		//init with -1
		pointBelongsTo[i] = -1;
	}
	
	//allocate memory for keeping track of sum of points that belong to each cluster
	double** pointsInCluster = (double **)malloc(numClusters*sizeof(double*));
	for(i=0;i<numClusters;i++) {
		pointsInCluster[i] = (double *)malloc(dim*sizeof(double));	//because we sum all the points for each cluster
	}
	//initialize #points count for each cluster
	int* numPointsInCluster = (int*)malloc(numClusters*sizeof(int));
	
	//ratio = (# points changing cluster)/(total # points)
	double ratio = threshold+1;
	int iter = 0;
	
	while(ratio > threshold && iter < maxIter) {
		for(i=0;i<numClusters;i++) {
			for(j=0;j<dim;j++) {
				pointsInCluster[i][j] = 0.0;
			}
			numPointsInCluster[i] = 0;
		}

		//for all points, calculate the closest cluster
		int totalPointsChanged = 0;
		for(i=0;i<numPoints;i++) {
			int closest_cluster = closest_cluster_calculator(points[i], cluster_centers, numClusters, dim);
			if(closest_cluster != pointBelongsTo[i]) {
				totalPointsChanged++;
				pointBelongsTo[i] = closest_cluster;
			}
			
			for(j=0;j<dim;j++) {
				pointsInCluster[closest_cluster][j] += points[i][j];
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
				cluster_centers[i][j] = pointsInCluster[i][j]/numPointsInCluster[i];
			}
		}
		
		iter++;
		ratio = ((double)totalPointsChanged)/numPoints;
	}
	
	printf("Iter=%d\n", iter);
	
	return cluster_centers;
}

int main(int argc, char* argv[]) {
	
	int numClusters = 5;
	double** init_centers = (double**)malloc(numClusters*sizeof(double*));
	int dim = 2;
	int totalNumPoints = 50000;
	
	
	///afs/andrew.cmu.edu/usr11/ndhruva/public/test.txt
	double** points = readFromFileForGP("/Users/neil/Documents/cluster.csv", dim, totalNumPoints);
	
	double** cluster_centers = get_cluster_centers_seq(points, totalNumPoints,
												  numClusters, dim,
												  1000000, totalNumPoints,
												  0.0);
//	int i, j;
//	//print cluster centers
//	for(i=0; i<numClusters; i++) {
//		for(j=0; j<dim; j++) {
//			printf("%f, ", cluster_centers[i][j]);
//		}
//		printf("\n");
//	}
	writeToFileForGP("./output/2doutput_seq.csv", cluster_centers, numClusters, dim);
	free(points);
	free(cluster_centers);
	return 0;
}
