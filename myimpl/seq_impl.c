#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include "allfunc.h"

/*
 * Calculates square of Euclidean distance between two points.
 */
double calc_squared_dist(double* point1, double* point2, int dim) {
	double retDist = 0.0f;
	int i=0;
	for(i=0; i<dim; i++) {
		retDist = retDist + pow((point2[i]-point1[i]), 2);
	}
	return retDist;
}

/*
 * Calculates the cluster center closest to a given points.
 */
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
				
	//create initial centers = first numClusters number of points
	double** cluster_centers = (double **)malloc(numClusters*sizeof(double*));
	int i, j;
	
	//code for assigning first numClusters points as initial centers
//	for(i=0;i<numClusters;i++) {
//		cluster_centers[i] = (double *)malloc(dim*sizeof(double));
//		for(j=0;j<dim;j++) {
//			cluster_centers[i][j] = points[i][j];			
//		}
//	}
	
	//assign random initial centers
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
	//initialize number of points count for each cluster
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
	
	printf("Total Iterations=%d\n", iter);
	
	return cluster_centers;
}

int main(int argc, char* argv[]) {
	
	int totalNumPoints = 0;
	int numClusters = 0;
	int dim = 2;						//default 2 dimensions
	double stopThreshold = 0.001;		//default threshold
	int maxIter = 100000;				//default max iterations
	char* inFile = "";
	char* outFile = "";
	int printTime = 0;
	extern char optopt;
	extern char* optarg;
	
	//process input
	char arg; int err;
	while((arg=getopt(argc,argv,"i:n:p:v::o::t::m::d::")) != -1) {
        switch (arg) {
            case 'i': {
				//input file
				inFile=optarg;
				break;
			}
			case 'n': {
				//number of clusters
				numClusters = atoi(optarg);
				break;
			}
			case 'p': {
				//total number of points
				totalNumPoints = atoi(optarg);
				break;
			}
			case 'o': {
				//output file
				outFile=optarg;
				break;
			}
            case 't': {
				//threshold for stopping (number of points changed/total number of points)
				stopThreshold=atof(optarg);
				break;
			}
			case 'm': {
				////maximum iterations (for stopping)
				maxIter=atoi(optarg);
				break;
			}
			case 'd': {
				//number of dimensions of data points
				dim=atoi(optarg);
				break;
			}
			case 'v': {
				//print the time taken for IO and Clustering (separately)
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
		return 0;
	}
	
	//malloc initial centers
	double** init_centers = (double**)malloc(numClusters*sizeof(double*));
	struct timeb tmb;
	
	ftime(&tmb);	
	double startInputTime = (double)tmb.time+(double)tmb.millitm/1000;
	//read input points from file
	double** points = readFromFileForGP(inFile, dim, totalNumPoints);
	ftime(&tmb);
	double endInputTime = (double)tmb.time+(double)tmb.millitm/1000;
	
	ftime(&tmb);
	double startClusteringTime = (double)tmb.time+(double)tmb.millitm/1000;
	//start clustering
	double** cluster_centers = get_cluster_centers_seq(points, totalNumPoints,
												  numClusters, dim,
												  maxIter, totalNumPoints,
												  stopThreshold);
	ftime(&tmb);			
	double endClusteringTime = (double)tmb.time+(double)tmb.millitm/1000;
	
	
	ftime(&tmb);
	double startOutputTime = (double)tmb.time+(double)tmb.millitm/1000;
	if(outFile && outFile[0] != '\0') {
		writeToFileForGP(outFile, cluster_centers, numClusters, dim);
	} else {
		printToTerminal(cluster_centers, numClusters, dim);
	}
	ftime(&tmb);
	double endOutputTime = (double)tmb.time+(double)tmb.millitm/1000;
	
	if(printTime == 1) {
		printf("IO Time: %lf\n", (endInputTime-startInputTime+endOutputTime-startOutputTime));
		printf("Clustering Time: %lf\n", (endClusteringTime-startClusteringTime));
		printf("Total Time: %lf\n", (endInputTime-startInputTime+endOutputTime-startOutputTime)+(endClusteringTime-startClusteringTime));
	}
	//free pointers
	free(points);
	free(cluster_centers);
	return 0;
}
