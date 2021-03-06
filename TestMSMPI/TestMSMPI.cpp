// TestMSMPI.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	//int provided = 0;
	//MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int length;
	char name[60];
	MPI_Get_processor_name(&name[0], &length);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		char helloStr[] = "Hello World!";
		MPI_Send(helloStr, _countof(helloStr), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		printf("Rank 0 sent string %s from rank 0\n", helloStr);
	}
	else if (rank == 1)
	{
		char helloStr[13];
		MPI_Recv(helloStr, _countof(helloStr), MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		printf("Rank 1 received string %s from rank 0\n", helloStr);
		MPI_Send(helloStr, _countof(helloStr), MPI_CHAR, 2, 0, MPI_COMM_WORLD);
	}
	else if (rank == 2)
	{
		char helloStr[13];
		MPI_Recv(helloStr, _countof(helloStr), MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		printf("Rank 2 received string %s from rank 1\n", helloStr);
	}

	MPI_Finalize();

    return 0;
}

