
#ifndef _allfunc_h
#define _allfunc_h

float calc_squared_dist(float*, float*, int);
int closest_cluster_calculator(float*, float**, int, int);
float** get_cluster_centers(float**, int, float**, int, int, int, float);

float** readFromFileForMPI(char*, int*, int, int, int, int, int, float**);
float** readFromFileForGP(char*, int, int);

#endif
