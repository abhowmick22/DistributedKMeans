/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         mpi_read.c                                                */
/*   File:         mpi_io.c                                                  */
/*   Description:  This program reads data points from a file and writes the */
/*                 clustering results to files using MPI-IO                  */
/*                 that implements a simple k-means clustering algorithm     */
/*   Input file format:                                                      */
/*                 ascii  file: each line contains 1 data object             */
/*                 binary file: first 4-byte integer is the number of data   */
/*                 objects and 2nd integer is the no. of features (or        */
/*                 coordinates) of each object                               */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*                                                                           */
/*   Copyright (C) 2005, Northwestern University                             */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
extern int errno;

#include <mpi.h>
#include "kmeans.h"


/*---< mpi_read() >----------------------------------------------------------*/
float** mpi_read(int       isBinaryFile,  /* flag: 0 or 1 */
                 char     *filename,      /* input file name */
                 int      *numObjs,       /* no. data objects (local) */
                 int      *numCoords,     /* no. coordinates */
                 MPI_Comm  comm)
{
    float    **objects;
    int        i, len, divd, rem;
    int        rank, nproc;
    MPI_Status status;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    
    else { /* ASCII format: let proc 0 read and distribute to others */
        if (rank == 0) {
            objects = file_read(0, filename, numObjs, numCoords);
            if (objects == NULL) *numObjs = -1;
        }

        /* broadcast global numObjs and numCoords to the rest proc */
        MPI_Bcast(numObjs,   1, MPI_INT, 0, comm);
        MPI_Bcast(numCoords, 1, MPI_INT, 0, comm);

        if (*numObjs == -1) {
            MPI_Finalize();
            exit(1);
        }
        if (*numObjs < nproc) {
            printf("Error: number of data points must be larger than the number of MPI processes.\n");
            MPI_Finalize();
            exit(1);
        }

        divd = (*numObjs) / nproc;
        rem  = (*numObjs) % nproc;

        if (rank == 0) {
            int index = (rem > 0) ? divd+1 : divd;

            /* index is the numObjs partitioned locally in proc 0 */
            (*numObjs) = index;

            /* distribute objects[] to other processes */
            for (i=1; i<nproc; i++) {
                int msg_size = (i < rem) ? (divd+1) : divd;
                MPI_Send(objects[index], msg_size*(*numCoords), MPI_FLOAT,
                         i, i, comm);
                index += msg_size;
            }

            /* reduce the objects[] to local size */
            objects[0] = realloc(objects[0],
                                 (*numObjs)*(*numCoords)*sizeof(float));
            assert(objects[0] != NULL);
            objects    = realloc(objects, (*numObjs)*sizeof(float*));
            assert(objects != NULL);
        }
        else {
            /*  local numObjs */
            (*numObjs) = (rank < rem) ? divd+1 : divd;

            /* allocate space for data points */
            objects    = (float**)malloc((*numObjs)            *sizeof(float*));
            assert(objects != NULL);
            objects[0] = (float*) malloc((*numObjs)*(*numCoords)*sizeof(float));
            assert(objects[0] != NULL);
            for (i=1; i<(*numObjs); i++)
                objects[i] = objects[i-1] + (*numCoords);

            MPI_Recv(objects[0], (*numObjs)*(*numCoords), MPI_FLOAT, 0,
                     rank, comm, &status);
        }
    }

    return objects;
}


/*---< mpi_write() >---------------------------------------------------------*/
int mpi_write(int        isOutFileBinary, /* flag: 0 or 1 */
              char      *filename,        /* input file name */
              int        numClusters,     /* no. clusters */
              int        numObjs,         /* no. data objects */
              int        numCoords,       /* no. coordinates (local) */
              float    **clusters,        /* [numClusters][numCoords] centers */
              int       *membership,      /* [numObjs] */
              int        totalNumObjs,    /* total no. data objects */
              MPI_Comm   comm,
              int        verbose)
{
    int        divd, rem, len, err;
    int        i, j, k, rank, nproc;
    char       outFileName[1024], fs_type[32], str[32], *delim;
    MPI_File   fh;
    MPI_Status status;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    delim = strchr(filename, ':');
    if (delim != NULL) {
        strncpy(fs_type, filename, delim-filename);
        fs_type[delim-filename] = '\0';
        /* real file name starts from delim+1 */
        delim++;
    }
    else
        delim = filename;

    /* output: the coordinates of the cluster centres ----------------------*/
    /* only proc 0 do this, because clusters[] are the same across all proc */
    if (rank == 0) {
        if (verbose) printf("Writing coordinates of K=%d cluster centers to file \"%s.cluster_centres\"\n",
                            numClusters, delim);
        sprintf(outFileName, "%s.cluster_centres", filename);
        err = MPI_File_open(MPI_COMM_SELF, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (err != MPI_SUCCESS) {
            char errstr[MPI_MAX_ERROR_STRING];
            int  errlen;
            MPI_Error_string(err, errstr, &errlen);
            printf("Error at opening file %s (%s)\n", outFileName,errstr);
            MPI_Finalize();
            exit(1);
        }
        if (isOutFileBinary) {
            MPI_File_write(fh, &numClusters, 1, MPI_INT, &status);
            MPI_File_write(fh, &numCoords,   1, MPI_INT, &status);
            MPI_File_write(fh, clusters[0], numClusters*numCoords, MPI_FLOAT, &status);
        }
        else {
            for (i=0; i<numClusters; i++) {
                sprintf(str, "%d ", i);
                MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                for (j=0; j<numCoords; j++) {
                    sprintf(str, "%f ", clusters[i][j]);
                    MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                }
                MPI_File_write(fh, "\n", 1, MPI_CHAR, &status);
            }
        }
        MPI_File_close(&fh);
    }

    /* output: the closest cluster centre to each of the data points --------*/
    if (rank == 0 && verbose)
        printf("Writing membership of N=%d data objects to file \"%s.membership\"\n",
               totalNumObjs, delim);

    if (isOutFileBinary) {
        sprintf(outFileName, "%s.membership", filename);
        err = MPI_File_open(comm, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (err != MPI_SUCCESS) {
            char errstr[MPI_MAX_ERROR_STRING];
            int  errlen;
            MPI_Error_string(err, errstr, &errlen);
            printf("Error at opening file %s (%s)\n", outFileName,errstr);
            MPI_Finalize();
            exit(1);
        }
        /* write numObjs from the 1st integer */
        if (rank == 0)
            MPI_File_write(fh, &totalNumObjs, 1, MPI_INT, &status);

        divd = totalNumObjs / nproc;
        rem  = totalNumObjs % nproc;
        len  = (rank < rem) ? rank*(divd+1) : rank*divd + rem;

        MPI_Offset disp = (len + 1) * sizeof(int);
        MPI_File_seek(fh, disp, MPI_SEEK_SET);
        MPI_File_write(fh, membership, numObjs, MPI_INT, &status);
        MPI_File_close(&fh);
    }
    else {
        if (rank == 0) { /* gather membership[] from all processes ----------*/
            int divd = totalNumObjs / nproc;
            int rem  = totalNumObjs % nproc;

            sprintf(outFileName, "%s.membership", filename);
            err = MPI_File_open(MPI_COMM_SELF, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
            if (err != MPI_SUCCESS) {
                char errstr[MPI_MAX_ERROR_STRING];
                int  errlen;
                MPI_Error_string(err, errstr, &errlen);
                printf("Error at opening file %s (%s)\n", outFileName,errstr);
                MPI_Finalize();
                exit(1);
            }

            /* first, print out local membership[] */
            for (j=0; j<numObjs; j++) {
                sprintf(str, "%d %d\n", j, membership[j]);
                MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
            }

            k = numObjs;
            for (i=1; i<nproc; i++) {
                numObjs = (i < rem) ? divd+1 : divd;
                MPI_Recv(membership, numObjs, MPI_INT, i, i, comm, &status);

                for (j=0; j<numObjs; j++) {
                    sprintf(str, "%d %d\n", k++, membership[j]);
                    MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                }
            }
            MPI_File_close(&fh);
        }
        else {
            MPI_Send(membership, numObjs, MPI_INT, 0, rank, comm);
        }
    }
    return 1;
}
