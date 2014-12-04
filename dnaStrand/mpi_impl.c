#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allfunc.h"
#include "mpi.h"

/* User defined type for weighted point */
typedef struct {
int a,t,g,c;		// weights for bases 'a', 't', 'g', 'c'
} weights;


/*
User defined function for reduce on dna strands
This function basically returns the mode of all elements
*/
void getWeights(weights *in, weights *inout, int *len,
		 MPI_Datatype *dptr)
{
	weights c;

	for (i=0; i< *len; ++i) {
		int j;
		c.a = in->a + inout->a;
		c.t = in->t + inout->t;
		c.g = in->g + inout->g;
		c.c = in->c + inout->c;
		*inout = c;
		in++; inout++;
	}	

}


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

char** get_cluster_centers_mpi(char** points, int numPoints,
							char** init_centers, int numClusters,
							int dim, int maxIter, 
							int totalNumPoints, double threshold) {

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
	
	//allocate memory for keeping track of base weights at each cluster
	weights** pointWeightsInCluster = (weights **)malloc(numClusters*sizeof(weights*));
	for(i=0;i<numClusters;i++) {
		pointWeightsInCluster[i] = (weights *)malloc(dim*sizeof(weights));	//because we accumulate weights all the points for each cluster
	}
	//initialize #points count for each cluster
	int* numPointsInCluster = (int*)malloc(numClusters*sizeof(int));
	
	//ratio = (# points changing cluster)/(total # points)
	double ratio = threshold+1;
	int iter = 0;	
	while(ratio > threshold && iter < maxIter) {
	
		for(i=0;i<numClusters;i++) {
			for(j=0;j<dim;j++) {
				// initialize weights
				pointWeightsInCluster[i][j].a = 0;			
				pointWeightsInCluster[i][j].t = 0;			
				pointWeightsInCluster[i][j].g = 0;			
				pointWeightsInCluster[i][j].c = 0;			
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
		
			/*	
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
			*/
			for(j=0;j<dim;j++){
				if(points[i][j] == 'a')				pointWeightsInCluster[closest_cluster][j].a += 1; 
				else if(points[i][j] == 't')		pointWeightsInCluster[closest_cluster][j].t += 1; 
				else if(points[i][j] == 'g')		pointWeightsInCluster[closest_cluster][j].g += 1; 
				else								pointWeightsInCluster[closest_cluster][j].c += 1; 
			}
			numPointsInCluster[closest_cluster]++;
		}

		int totalPointsChanged = 0;
		//int totalNumPoints = 0;
		int* tempNumPointsInCluster = (int*)malloc(numClusters*sizeof(int));
		
		//get the total sum of number of points that changed centers
		MPI_Allreduce(&localPointsChanged, &totalPointsChanged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//get the total number of points in the dataset (could also be passed as a param to this function)
		//MPI_Allreduce(&numPoints, &totalNumPoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//for each center, get total number of points belonging to it
		for(i=0;i<numClusters;i++) {
			MPI_Allreduce(&numPointsInCluster[i], &tempNumPointsInCluster[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		}
		int* freeNumPointsInCluster = numPointsInCluster;
		numPointsInCluster = tempNumPointsInCluster;
		free(freeNumPointsInCluster);

		//to compute each cluster center, get the weights of all coordinates of all points belonging to it
		weights** tempPointWeightsInCluster = (weights **)malloc(numClusters*sizeof(weights*));
		temp = (weights*)malloc(numClusters*dim*sizeof(weights));
		for(i=0;i<numClusters;i++) {
			tempPointWeightsInCluster[i] = &temp[i*dim];
		}

		// create the getWeights user-op
		MPI_Op myOp;
		MPI_Op_create(getWeights, True, &myOp);

		// explain to MPI how type weights is defined
		MPI_Datatype wtype;
        MPI_Type_contiguous(4, MPI_INT, &wtype);
        MPI_Type_commit(&wtype);

		for(i=0;i<numClusters;i++) {
			// TODO: How to do this step ?, use mpi_op_create
			MPI_Allreduce(pointWeightsInCluster[i], tempPointWeightsInCluster[i], dim, wtype, myOp, MPI_COMM_WORLD);
		}
		char** freePointWeightsInCluster = pointWeightsInCluster;
		pointWeightsInCluster = tempPointWeightsInCluster;
		free(freePointWeightsInCluster);
		
		//assign new cluster centers
		for(i=0; i<numClusters; i++) {
			if(numPointsInCluster[i]==0) {
				//shouldn't happen because we start with centers chosen from input data points
				continue;
			}
			for(j=0; j<dim; j++) {
				//cluster_centers[i][j] = pointsInCluster[i][j];		// the base at this position is already majority, if it exists
				// find the most frequent base	
				int maxIndex = 0;
				int max = 0;
				if(pointWeightsInCluster[i][j].a > max){
					max = pointWeightsInCluster[i][j].a;
					maxIndex=0;
				}
				if(pointWeightsInCluster[i][j].t > max){
					max = pointWeightsInCluster[i][j].t;
					maxIndex=1;
				}
				if(pointWeightsInCluster[i][j].g > max){
					max = pointWeightsInCluster[i][j].g;
					maxIndex=2;
				}
				if(pointWeightsInCluster[i][j].c > max){
					max = pointWeightsInCluster[i][j].c;
					maxIndex=3;
				}
				
				// set the centroid
				if(maxIndex == 0)		cluster_centers[i][j] = 'a';
				else if(maxIndex == 1)	cluster_centers[i][j] = 't';
				else if(maxIndex == 2)	cluster_centers[i][j] = 'g';
				else					cluster_centers[i][j] = 'c';
			}
		}
		
		iter++;
		ratio = ((double)totalPointsChanged)/totalNumPoints;
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
	
	int numClusters = 5;
	int* numPoints = (int*)malloc(sizeof(int));
	char** init_centers = (char**)malloc(numClusters*sizeof(char*));
	int dim = 2;
	int totalNumPoints = 50000;
	
	///afs/andrew.cmu.edu/usr11/ndhruva/public/test.txt
	char** points = readFromFileForMPI("/home/abhishek/15-640/project4/dnaStrand/cluster.csv", numPoints,
					   					dim, totalNumPoints, numClusters,
					   					rank, numProcs, init_centers);
				
	MPI_Barrier(MPI_COMM_WORLD);
	double endInputTime = MPI_Wtime();
				
	if(points==NULL) {
		MPI_Finalize();
		return 0;
	}
		
	double startClusteringTime = MPI_Wtime();
	char** cluster_centers = get_cluster_centers_mpi(points, *numPoints,
												  init_centers, numClusters,
												  dim, 1000000, 
												  totalNumPoints, 0.0);
	MPI_Barrier(MPI_COMM_WORLD);
	double endClusteringTime = MPI_Wtime();
	
	if(rank==0) {
		double startOutputTime = MPI_Wtime();
		writeToFileForGP("./output/dnaOutput_mpi.csv", cluster_centers, numClusters, dim);
		double endOutputTime = MPI_Wtime();
		
		//total IO time
		printf("Total IO Time = %f\n", (endInputTime-startInputTime)+(endOutputTime-startOutputTime));
		//total Clustering time
		printf("Total Clustering Time = %f\n", (endClusteringTime-startClusteringTime));
	}
	
	
	free(points);
	free(cluster_centers);
	free(init_centers);
	free(numPoints);
	
	MPI_Finalize();
	return 0;
}
