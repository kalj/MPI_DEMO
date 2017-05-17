/* -*- c-basic-offset:4; tab-width:4; indent-tabs-mode:nil -*-
 *
 * @(#)Demo_pio_mat.c
 * @author Karl Ljungkvist <k.ljungkvist@gmail.com>
 * Date: 2017-05-17
 */

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    int np, rank;
    int p;
    int locsize;
    int size = 16;
    int coord[2];
    int dims[2];
    int matdims[2];
    int locmatdims[2];
    int starts[2];
    MPI_Datatype subarray;
    MPI_File fh;
    MPI_Status status;

    MPI_Init(&argc,&argv);

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) == NULL){
        perror("getcwd() error");
        MPI_Finalize();
        return 0;
    }
    char full_path[1024];
    sprintf (full_path, "%s/mat.dat", cwd);
    const char * file =  (const char * )  &full_path[0];

    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    p = round(sqrtf(np));

    locsize = size/p;
    coord[0] = rank/p;
    coord[1] = rank%p;
    dims[0] = p;
    dims[1] = p;
    matdims[0] = size;
    matdims[1] = size;
    locmatdims[0] = locsize;
    locmatdims[1] = locsize;

    starts[0] = coord[0]*locsize;
    starts[1] = coord[1]*locsize;

    char *local_array = malloc(locsize*locsize*sizeof(char));

    for(int i = 0; i < locsize; ++i) {
        for(int j = 0; j < locsize; ++j) {
            local_array[i*locsize + j] = 'a'+rank;
        }
    }

    MPI_Type_create_subarray(2,matdims,locmatdims, starts, MPI_ORDER_C,MPI_CHAR,&subarray);
    MPI_Type_commit(&subarray);

    MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_WRONLY | MPI_MODE_CREATE , MPI_INFO_NULL, &fh );

    MPI_File_set_view(fh, 0, MPI_CHAR, subarray, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, local_array, locsize*locsize, MPI_CHAR, &status);
    MPI_File_close(&fh);

    MPI_Finalize();

    return 0;
}
/*
-bash-4.1$ mpicc -o Demopiomat ../Demo_pio_mat.c -lm
-bash-4.1$ mpirun -n 4 Demopiomat
-bash-4.1$ fmt mat.dat | fold -w 16
aaaaaaaabbbbbbbb
aaaaaaaabbbbbbbb
aaaaaaaabbbbbbbb
aaaaaaaabbbbbbbb
aaaaaaaabbbbbbbb
aaaaaaaabbbbbbbb
aaaaaaaabbbbbbbb
aaaaaaaabbbbbbbb
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
ccccccccdddddddd
*/
