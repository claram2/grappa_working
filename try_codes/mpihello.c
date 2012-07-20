#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(){
MPI_Init(NULL, NULL);
int comm_sz;
int my_rank;
int i;
char greeting[100];
MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

sprintf(greeting, "hello from process %d of %d\n", my_rank, comm_sz);

MPI_Send(greeting, 100, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

if(my_rank == 0){
  for(i = 0; i < comm_sz; ++i){
    MPI_Recv(greeting, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf(greeting);
  }
}
MPI_Finalize();
}




