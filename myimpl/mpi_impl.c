#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allfunc.h"
#include "mpi.h"

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

double** get_cluster_centers_mpi(double** points, int numPoints,
								double** init_centers, int numClusters,
								int dim, int maxIter,
								int totalNumPoints, double threshold) {

	int i, j;
	double* temp;
	//initialize cluster centers with init_centers
	double** cluster_centers = (double **)malloc(numClusters*sizeof(double*));
	temp = (double*)malloc(numClusters*dim*sizeof(double));
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
	
	//allocate memory for keeping track of sum of points that belong to each cluster
	double** pointsInCluster = (double **)malloc(numClusters*sizeof(double*));
	temp = (double*)malloc(numClusters*dim*sizeof(double));
	for(i=0;i<numClusters;i++) {
		pointsInCluster[i] = &temp[i*dim];	//because we sum all the points for each cluster;
										    //also, contiguous memory very important for allreduce
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
		int localPointsChanged = 0;
		for(i=0;i<numPoints;i++) {
			int closest_cluster = closest_cluster_calculator(points[i], cluster_centers, numClusters, dim);
			if(closest_cluster != pointBelongsTo[i]) {
				localPointsChanged++;
				pointBelongsTo[i] = closest_cluster;
			}			
			for(j=0;j<dim;j++) {
				pointsInCluster[closest_cluster][j] += points[i][j];
			}
			numPointsInCluster[closest_cluster]++;
		}
		
		int totalPointsChanged = 0;
		//int totalNumPoints = 0;
		int* tempNumPointsInCluster = (int*)malloc(numClusters*sizeof(int));
		
		//get the total sum of number of points that changed centers
		MPI_Allreduce(&localPointsChanged, &totalPointsChanged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//for each center, get total number of points belonging to it
		MPI_Allreduce(numPointsInCluster, tempNumPointsInCluster, numClusters, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		int* freeNumPointsInCluster = numPointsInCluster;
		numPointsInCluster = tempNumPointsInCluster;
		free(freeNumPointsInCluster);
		//for each center, get the sum of all coordinates of all points belonging to it
		double** tempPointsInCluster = (double **)malloc(numClusters*sizeof(double*));
		temp = (double*)malloc(numClusters*dim*sizeof(double));
		for(i=0;i<numClusters;i++) {
			tempPointsInCluster[i] = &temp[i*dim];
		}
		MPI_Allreduce(pointsInCluster[0], tempPointsInCluster[0], numClusters*dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		double** freePointsInCluster = pointsInCluster;
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
		ratio = ((double)totalPointsChanged)/totalNumPoints;
	}
	
	return cluster_centers;
	
}

int main(int argc, char* argv[]) {
	
	int totalNumPoints = 0;
	int numClusters = 0;
	int dim = 2;						//default 2 dimensions
	char* inFile = "";				//"/Users/neil/Documents/cluster.csv", /afs/andrew.cmu.edu/usr11/ndhruva/public/test.txt
	char* outFile = "";				//"./output/2doutput_mpi.csv"
	double stopThreshold = 0.001;		//default threshold
	int maxIter = 100000;				//default max iterations
	int printTime = 0;
	extern char optopt;
	extern char* optarg;
	
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
	
	
	//process input
	char arg; int err;
	while((arg=getopt(argc,argv,"i:n:p:otmdv")) != -1) {
        switch (arg) {
            case 'i': {
				inFile=optarg;
				break;
			}
			case 'n': {
				numClusters = atoi(optarg);
				break;
			}
			case 'p': {
				totalNumPoints = atoi(optarg);
				break;
			}
			case 'o': {
				outFile=optarg;
				break;
			}
            case 't': {
				stopThreshold=atof(optarg);
				break;
			}
			case 'm': {
				maxIter=atoi(optarg);
				break;
			}
			case 'd': {
				dim=atoi(optarg);
				break;
			}
			case 'v': {
				printTime=atoi(optarg);
				break;
			}
			case '?': {
				fprintf(stderr, "Unrecognized option -%c.\n", optopt);
				err = 1;
				break;
			}
            default: {
				break;
			}
        }
    }
	
	if(err==1 || (inFile && inFile[0] == '\0') || numClusters <= 0 || totalNumPoints <= 0 || stopThreshold <= 0.0 || maxIter <= 0) {
		MPI_Finalize();
		return 0;
	}
	
	int* numPoints = (int*)malloc(sizeof(int));
	double** init_centers = (double**)malloc(numClusters*sizeof(double*));
		
	//block for all to synchronize
	MPI_Barrier(MPI_COMM_WORLD);
	double startInputTime = MPI_Wtime();
	
	
	double** points = readFromFileForMPI(inFile, numPoints,
					   					dim, totalNumPoints, numClusters,
					   					rank, numProcs, init_centers);
	
	MPI_Barrier(MPI_COMM_WORLD);
	double endInputTime = MPI_Wtime();
				
	if(points==NULL) {
		MPI_Finalize();
		return 0;
	}
		
	double startClusteringTime = MPI_Wtime();
	double** cluster_centers = get_cluster_centers_mpi(points, *numPoints,
												  init_centers, numClusters,
												  dim, maxIter,
												  totalNumPoints, stopThreshold);
	
	MPI_Barrier(MPI_COMM_WORLD);
	double endClusteringTime = MPI_Wtime();
	
	if(rank==0) {
		double startOutputTime = MPI_Wtime();
		if(outFile && outFile[0] != '\0') {
			writeToFileForGP(outFile, cluster_centers, numClusters, dim);
		} else {
			printToTerminal(cluster_centers, numClusters, dim);
		}
		double endOutputTime = MPI_Wtime();
		if(printTime == 1) {
			//total IO time
			printf("Total IO Time = %f\n", (endInputTime-startInputTime)+(endOutputTime-startOutputTime));
			//total Clustering time
			printf("Total Clustering Time = %f\n", (endClusteringTime-startClusteringTime));
		}
	}
	
	
	free(points);
	free(cluster_centers);
	free(init_centers);
	free(numPoints);
	
	MPI_Finalize();
	return 0;
}

