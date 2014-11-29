//
//  allfunc.h
//  
//
//  Created by Neil on 11/28/14.
//
//

#ifndef _allfunc_h
#define _allfunc_h

float calc_squared_dist(float*, float*, int);
int closest_cluster_calculator(float*, float**, int, int);
float** get_cluster_centers(float** points, int numPoints, int numClusters, int dim, int maxIter, float threshold);


#endif
