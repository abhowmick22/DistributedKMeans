
#ifndef _allfunc_h
#define _allfunc_h

int calc_string_distance(char*, char*, int);
int closest_cluster_calculator(char*, char**, int, int);
char** get_cluster_centers(char**, int, char**, int, int, int, float);

char** readFromFileForMPI(char*, int*, int, int, int, int, int, char**);
char** readFromFileForGP(char*, int, int);

#endif
