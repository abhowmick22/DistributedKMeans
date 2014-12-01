#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allfunc.h"
#include "mpi.h"

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
							char** init_centers, int numClusters,
							int dim, int maxIter, float threshold) {

	int i, j;
	char* temp;
	//initialize cluster centers with init_centers
	char** cluster_centers = (char **)malloc(numClusters*sizeof(char*));
	temp = (char*)malloc(numClusters*dim*sizeof(char));
	for(i=0;i<numClusters;i++) {
		cluster_centers[i] = &temp[i*dim];
	}
	for(i=0;i<numClusters;i++) {
		for(j=0;j<dim;j++) {
			cluster_centers[i][j] = init_centers[i][j];
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
		
		
		//TODO : Beyond this


		int totalPointsChanged = 0;
		int totalNumPoints = 0;
		int* tempNumPointsInCluster = (int*)malloc(numClusters*sizeof(int));
		
		//get the total sum of number of points that changed centers
		MPI_Allreduce(&localPointsChanged, &totalPointsChanged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//get the total number of points in the dataset (could also be passed as a param to this function)
		MPI_Allreduce(&numPoints, &totalNumPoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//for each center, get total number of points belonging to it
		for(i=0;i<numClusters;i++) {
			MPI_Allreduce(&numPointsInCluster[i], &tempNumPointsInCluster[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		}
		int* freeNumPointsInCluster = numPointsInCluster;
		numPointsInCluster = tempNumPointsInCluster;
		free(freeNumPointsInCluster);
		//for each center, get the sum of all coordinates of all points belonging to it
		float** tempPointsInCluster = (float **)malloc(numClusters*sizeof(float*));
		temp = (float*)malloc(numClusters*dim*sizeof(float));
		for(i=0;i<numClusters;i++) {
			tempPointsInCluster[i] = &temp[i*dim];
		}
		for(i=0;i<numClusters;i++) {
			MPI_Allreduce(pointsInCluster[i], tempPointsInCluster[i], dim, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
		}
		float** freePointsInCluster = pointsInCluster;
		pointsInCluster = tempPointsInCluster;
		free(freePointsInCluster);
		
		//assign new cluster centers
		for(i=0; i<numClusters; i++) {
			if(numPointsInCluster[i]==0) {
				//shouldn't happen because we start with centers chosen from input data points
				continue;
			}
			for(j=0; j<dim; j++) {
				cluster_centers[i][j] = pointsInCluster[i][j]/numPointsInCluster[i];
			}
		}
		
		iter++;
		ratio = ((float)totalPointsChanged)/totalNumPoints;
	}
	
	return cluster_centers;
	
}

int main(int argc, char* argv[]) {
	
	int flag;
	MPI_Initialized(&flag);

	int numProcs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	
	if(!flag) {
		MPI_Init(&argc, &argv);
	}
	
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);
	
	int numClusters = 3;
	int* numPoints = (int*)malloc(sizeof(int));
	float** init_centers = (float**)malloc(numClusters*sizeof(float*));
	int dim = 2;
	int totalNumPoints = 9;
	
	///afs/andrew.cmu.edu/usr11/ndhruva/public/test.txt
	float** points = readFromFileForMPI("/home/abhishek/15-640/project4/dnaStrand/test.txt", numPoints,
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
		//print cluster centers
		for(i=0; i<numClusters; i++) {
			for(j=0; j<dim; j++) {
				printf("%f, ", cluster_centers[i][j]);
			}
			printf("\n");
		}
	}
	free(points);
	free(cluster_centers);
	free(init_centers);
	free(numPoints);
	
	MPI_Finalize();
	return 0;
}
