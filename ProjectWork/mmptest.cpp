#include <stdio.h>
#include "mpi.h"
#include <math.h>


int main(int argc,char* argv[]) {

	int myid, numprocs;
	MPI_Init(&argc,&argv);
	printf("Loading the Data...\n");
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	printf("my id is: %d\n",myid);
			



	MPI_Finalize();
return 0;
}
