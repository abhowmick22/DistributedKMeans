#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include "allfunc.h"

/* User defined type for weighted point */
typedef struct {
int a,t,g,c;		// weights for bases 'a', 't', 'g', 'c'
} weights;


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

char** get_cluster_centers_seq(char** points, int numPoints,
							int numClusters, int dim,
							int maxIter, int totalNumPoints, 
							double threshold) {
				
	//create initial centers = first numClusters number of
	//points
	char** cluster_centers = (char **)malloc(numClusters*sizeof(char*));
	int i, j;
/*
	for(i=0;i<numClusters;i++) {
		cluster_centers[i] = (char *)malloc(dim*sizeof(char));
		for(j=0;j<dim;j++) {
			cluster_centers[i][j] = points[i][j];			
		}
	}
*/

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
		cluster_centers[i] = (char *)malloc(dim*sizeof(char));
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
		int totalPointsChanged = 0;
		for(i=0;i<numPoints;i++) {
			int closest_cluster = closest_cluster_calculator(points[i], cluster_centers, numClusters, dim);
			if(closest_cluster != pointBelongsTo[i]) {
				totalPointsChanged++;
				pointBelongsTo[i] = closest_cluster;
			}
			
			// Majority voting algorithm for all positions, for each incoming data point 
			// CHECK: counter never goes below zero
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
		
		//assign new cluster centers
		for(i=0; i<numClusters; i++) {
			if(numPointsInCluster[i]==0) {
				//shouldn't happen
				continue;
			}
			for(j=0; j<dim; j++) {
				//cluster_centers[i][j] = pointsInCluster[i][j];	// the base at this position is already majority, if it exists
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
		ratio = ((double)totalPointsChanged)/numPoints;
	}
	
	printf("Iter=%d\n", iter);
	
	return cluster_centers;
}

int main(int argc, char* argv[]) {
	
	int numClusters = 5;
	char** init_centers = (char**)malloc(numClusters*sizeof(char*));
	int dim = 2;
	int totalNumPoints = 50000;
	struct timeb tmb;
	
	ftime(&tmb);	
	double startInputTime = (double)tmb.time+(double)tmb.millitm/1000;
	///afs/andrew.cmu.edu/usr11/ndhruva/public/test.txt
	char** points = readFromFileForMPI("/home/abhishek/15-640/project4/dnaStrand/cluster.csv", dim, totalNumPoints);
	ftime(&tmb);
	double endInputTime = (double)tmb.time+(double)tmb.millitm/1000;
	printf("Read done\n");
	
	ftime(&tmb);
	double startClusteringTime = (double)tmb.time+(double)tmb.millitm/1000;
	char** cluster_centers = get_cluster_centers_seq(points, totalNumPoints,
												  numClusters, dim, 
												  1000000, totalNumPoints, 
												  0.0);
	ftime(&tmb);			
	double endClusteringTime = (double)tmb.time+(double)tmb.millitm/1000;
	printf("Clustering done\n");
	
	
	ftime(&tmb);
	double startOutputTime = (double)tmb.time+(double)tmb.millitm/1000;
	writeToFileForGP("./output/dnaOutput_seq.csv", cluster_centers, numClusters, dim);
	ftime(&tmb);
	double endOutputTime = (double)tmb.time+(double)tmb.millitm/1000;
	printf("Write done\n");
	
	
	printf("IO time: %lf\n", (endInputTime-startInputTime+endOutputTime-startOutputTime));
	printf("Clustering time: %lf\n", (endClusteringTime-startClusteringTime));
	
	free(points);
	free(cluster_centers);
	return 0;
}
