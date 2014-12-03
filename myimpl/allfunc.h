
#ifndef _allfunc_h
#define _allfunc_h

double calc_squared_dist(double*, double*, int);
int closest_cluster_calculator(double*, double**, int, int);
double** get_cluster_centers_mpi(double**, int, double**, int, int, int, int, double);
double** get_cluster_centers_seq(double**, int, int, int, int, int, double);

double** readFromFileForMPI(char*, int*, int, int, int, int, int, double**);
double** readFromFileForGP(char*, int, int);
void writeToFileForGP(char*, double**, int, int);

void plotResults(double*, double*, int);

#endif
