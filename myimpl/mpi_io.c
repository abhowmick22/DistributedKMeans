#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "allfunc.h"
#include "mpi.h"

/*
 * Return all points after reading from file and sending 
 * points to their corresponding processors for nearest cluster
 * center calculation.
 */
float** readFromFileForMPI(char* filename, int* numPoints,
							int dim, int totalNumPoints, int numClusters,
							int rank, int numProcs, float** init_centers) {
				
	MPI_Status status;
	float* temp;
	float** points;
	int sigCont;
	//logic: the master reads the file using general_io and sends
	//certain number of points to each processor to calculate
	//the nearest cluster center of each point
	
	if(totalNumPoints < numProcs) {
		//can't have this
		printf("Cannot have more processors than number of input points. Exiting..\n");
		//quit
		MPI_Finalize();
		exit(1);
	}
	
	//master reads the points from file using general_io
	if(rank==0) {
		points = readFromFileForGP(filename, dim, totalNumPoints);
		if(points == NULL) {
			//signal all processes to stop
			sigCont = 0;
		} else {
			//signal all processes to continue
			sigCont = 1;
		}
		
	}
	MPI_Bcast(&sigCont, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(!sigCont) {
		//quit
		MPI_Finalize();
		exit(1);
	}
	
	//now the master can distribute the points read to different processes
	//logic: distribute totalNumPoints/numProcs points to all.
	//but, send 1 more point to certain procs when totalNumPoints%numProcs!=0
	//This equally distributes points.
	int pointsToAll = totalNumPoints/numProcs;
	int pointToSome = totalNumPoints%numProcs;	//extra point is sent to those procs with
												//rank < remaining #points (rank starts at 0)
	int i, j;
	temp = (float*)malloc(numClusters*dim*sizeof(float));
	for(i=0; i<numClusters; i++) {
		init_centers[i] = &temp[i*dim];
	}
	
	if(rank == 0) {
		//master node
		
		//select initial centers randomly
		int* prevCenters = (int*)malloc(numClusters*sizeof(int));		
		for(i=0; i<numClusters; i++) {
			prevCenters[i] = -1;
		}
		for(i=0;i<numClusters;) {
			//srand(time(NULL));	//TODO: can be remove to reduce delay
			int r = i;//TODO: rand() % totalNumPoints;
			//check if the same has been encountered before
			for(j=0;j<i;j++) {
				if(prevCenters[j]==r) {
					break;
				}
			}
			if(j!=i) {
				continue;
			}
			prevCenters[i] = r;
			for(j=0;j<dim;j++) {
				init_centers[i][j] = points[r][j];
			}
			i++;
		}
		free(prevCenters);
		
		//broadcast init centers to all other nodes
		MPI_Bcast(init_centers[0], numClusters*dim, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		int startIndex = 0;
		//points to be analysed by master
		if(pointToSome != 0) {
			startIndex = pointsToAll+1;
		} else {
			startIndex = pointsToAll;
		}
		
		(*numPoints) = startIndex;
		
		//distribute the remaining to other procs
		for(j=1;j<numProcs;j++) {
			int numPointsToSend = pointsToAll;
			if(j < pointToSome) {
				numPointsToSend++;
			}
			//tag each message with the rank of the receiver
			MPI_Send(points[startIndex], numPointsToSend*dim, MPI_FLOAT, j, j, MPI_COMM_WORLD);
			startIndex += numPointsToSend;
		}
		
		//now, the master should process only the first part of the points
		float** newPoints = (float**)malloc((*numPoints)*sizeof(float*));
		temp = malloc((*numPoints)*dim*sizeof(float));
		for(i=0; i<(*numPoints); i++) {
			newPoints[i] = &temp[i*dim];
		}
		for(i=0;i<(*numPoints);i++) {
			for(j=0;j<dim;j++) {
				newPoints[i][j] = points[i][j];
			}
		}
		float** oldPoints = points;
		points = newPoints;
		free(oldPoints);
	} else {
		//helper node		
		
		(*numPoints) = pointsToAll;
		if(rank < pointToSome) {
			(*numPoints)++;
		}
		
		points = (float**)malloc((*numPoints)*sizeof(float*));
		temp = (float*)malloc((*numPoints)*dim*sizeof(float));
		for(i=0; i<(*numPoints); i++) {
			points[i] = &temp[i*dim];
		}
		MPI_Recv(points[0], (*numPoints)*dim, MPI_FLOAT, 0, rank, MPI_COMM_WORLD, &status);
	}
	
	
	return points;
	
}
