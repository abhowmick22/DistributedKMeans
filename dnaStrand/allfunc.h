
#ifndef _allfunc_h
#define _allfunc_h

int calc_string_distance(char*, char*, int);
int closest_cluster_calculator(char*, char**, int, int);
char** get_cluster_centers_mpi(char**, int, char**, int, int, int, int, double);
char** get_cluster_centers_seq(char**, int, int, int, int, int, double);

char** readFromFileForMPI(char*, int*, int, int, int, int, int, char**);
char** readFromFileForGP(char*, int, int);
void writeToFileForGP(char*, char**, int, int);
void printToTerminal(char**, int, int);

#endif
