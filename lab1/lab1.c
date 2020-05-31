#include <stdio.h>
#include <mpi.h>

int main(int argc,char *argv[]) 
{
    int i;
    int rank;
    int size;
    int j = 0;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(id != 0) {
        MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } else if(rank == 0) {
        while(j<p-1) {
            MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Hello rank %d.I'm rank %d \n",rank,i);
            fflush(stdin);
            j++;
        }
    }
    MPI_Finalize();
    
    return 0;
}
