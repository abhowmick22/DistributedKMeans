#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allfunc.h"

float calc_squared_dist(float* point1, float* point2, int dim) {
	float retDist = 0.0f;
	int i=0;
	for(i=0; i<dim; i++) {
		retDist = retDist + pow((point2[i]-point1[i]), 2);
	}
	return retDist;
}

int closest_cluster_calculator(float* point, float** cluster_centers,
								int numClusters, int dim) {
	float min_dist = calc_squared_dist(point, cluster_centers[0], dim);
	int closest_center = 0;
	int i=1;
	for(; i<numClusters; i++) {
		float temp_dist = calc_squared_dist(point, cluster_centers[i], dim);
		if(temp_dist < min_dist) {
			min_dist = temp_dist;
			closest_center = i;
		}
	}
	return closest_center;
}

float** get_cluster_centers(float** points, int numPoints,
							int numClusters, int dim,
							int maxIter, float threshold) {
				
	//create initial centers = first numClusters number of
	//points
	float** cluster_centers = (float **)malloc(numClusters*sizeof(float*));
	int i, j;
	for(i=0;i<numClusters;i++) {
		cluster_centers[i] = (float *)malloc(dim*sizeof(float));
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
	
	//allocate memory for keeping track of sum of points that belong to each cluster
	float** pointsInCluster = (float **)malloc(numClusters*sizeof(float*));
	for(i=0;i<numClusters;i++) {
		pointsInCluster[i] = (float *)malloc(dim*sizeof(float));	//because we sum all the points for each cluster
	}
	//initialize #points count for each cluster
	int* numPointsInCluster = (int*)malloc(numClusters*sizeof(int));
	
	//ratio = (# points changing cluster)/(total # points)
	float ratio = threshold+1;
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
		ratio = ((float)totalPointsChanged)/numPoints;
	}
	
	printf("Iter=%d\n", iter);
	
	return cluster_centers;
}

int main(int argc, char* argv[]) {
	
	int numClusters = 3;
	int* numPoints = (int*)malloc(sizeof(int));
	float** init_centers = (float**)malloc(numClusters*sizeof(float*));
	int dim = 6;
	int totalNumPoints = 18;
	
	
	///afs/andrew.cmu.edu/usr11/ndhruva/public/test.txt
	float** points = readFromFileForMPI("/Users/neil/Documents/test.txt", numPoints,
					   					dim, totalNumPoints, numClusters,
					   					rank, numProcs, init_centers);
	
	if(points==NULL) {
		MPI_Finalize();
		return 0;
	}
	float** cluster_centers = get_cluster_centers(points, *numPoints,
												  init_centers, numClusters,
												  dim, 1000000, 0.0);
	int i, j;
	if(rank==0) {
		printf("\n");
		for(i=0; i<numClusters; i++) {
			for(j=0; j<dim; j++) {
				printf("%f, ", init_centers[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		//print cluster centers
		for(i=0; i<numClusters; i++) {
			for(j=0; j<dim; j++) {
				printf("%f, ", cluster_centers[i][j]);
			}
			printf("\n");
		}
	}
	/*
	 int i, j;
	 printf("Rank: %d, NumPoints: %d\n", rank, *numPoints);
	 for(i=0;i<*numPoints;i++) {
	 for(j=0;j<dim;j++) {
	 printf("%f, ", points[i][j]);
	 }
	 printf("\n");
	 }
	 */
	free(points);
	free(cluster_centers);
	free(init_centers);
	free(numPoints);
//	
//	float** points = (float **) malloc(3*sizeof(float*));
//	points[0] = (float *) malloc(3*sizeof(float));
//	points[1] = (float *) malloc(3*sizeof(float));
//	points[2] = (float *) malloc(3*sizeof(float));
//	points[0][0]=2.0f;
//	points[0][1]=2.0f;
//	points[0][2]=2.0f;
//	points[1][0]=1.0f;
//	points[1][1]=1.0f;
//	points[1][2]=1.0f;
//	points[2][0]=1.5f;
//	points[2][1]=1.5f;
//	points[2][2]=1.5f;
//	//printf("%f\n", calc_squared_dist(point1, point2, 3));
//	float** cluster_centers = get_cluster_centers(points, 3, 3, 3, 100, 0.01);
//	int i=0, j=0;
//	for(i=0;i<3;i++) {
//		for(j=0; j<3; j++) {
//			printf("%f, ", cluster_centers[i][j]);
//		}
//		printf("\n");
//	}
//				
//				
//	free(points);
	return 0;
}
