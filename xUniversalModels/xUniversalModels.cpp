// xUniversalModels.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "ComponentCollection.h"
#include "Flash.h"
#include "Util.h"
#include "ODESolver.h"

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	realtype theta = NV_Ith_S(y, 0);
	realtype omega = NV_Ith_S(y, 1);
	realtype omegap = -sin(theta);
	NV_Ith_S(ydot, 0) = omega;
	NV_Ith_S(ydot, 1) = omegap;
	return 0;
}

extern "C"
{
	__declspec(dllexport) void RunUnivModel(TemperatureProfileEnum reactorTemperatureProfile, int ysize, double* yin, double* yout, int tsize, double* temperatureProfile)
	{
		ODESolver* ode = new ODESolver();
		ode->SetReactorData();

		int n = sizeof(yin);
		bool isOK = ode->RunReactorModel(reactorTemperatureProfile, ysize, yin, tsize, temperatureProfile);
	}

	__declspec(dllexport) void TestMPI(int size, int *a, int *b)
	{
		ODESolver* ode = new ODESolver();
		ode->RunA();

		double A = 10.0;
		double B = -10.0;
		double Tol = 0.00001;
		double x = Util::brents_fun(Util::TestFunction, A, B, Tol, 1000);
		x = Util::TestFunction(x);

		MPI_Init(NULL, NULL);

		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank == 0)
		{
			char helloStr[] = "Hello World!";
			//MPI_Send(helloStr, _countof(helloStr), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
			MPI_Send(helloStr, 13, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
			//printf("Rank 0 sent string %s from rank 0\n", helloStr);
		}
		else if (rank == 1)
		{
			char helloStr[13];
			MPI_Recv(helloStr, _countof(helloStr), MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			//printf("Rank 1 received string %s from rank 0\n", helloStr);
		}

		MPI_Finalize();
	}

	__declspec(dllexport) int add(int a, int b)
	{
		return a + b;
	}

	__declspec(dllexport) int subtract(int a, int b)
	{
		return a - b;
	}

	__declspec(dllexport) void test1(int a[], int b[])
	{
		int size = sizeof(a);
		for (int i = 0; i < size; ++i)
			a[i] = a[i] - b[i];
	}

	__declspec(dllexport) bool test2(int a[2][3])
	{
		int rows = sizeof a / sizeof a[0]; // 2 rows  
		int cols = sizeof a[0] / sizeof(int); // 5 cols
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				a[i][j] += 1;
			}
		}

		return false;
	}

	__declspec(dllexport) void test3(int size, double *a, double *b)
	{
		int n = sizeof(a);
		n = sizeof(double);

		for (int i = 0; i < size; ++i)
			*(a + i) = *(a + i) - *(b + i);
	}

	__declspec(dllexport) void test4(UserDataStruct* ud)
	{
		ud->n += 1;
		ud->x += 1.0;
	}

	__declspec(dllexport) void test5(char* strx, _Out_ LPTSTR lpString, _In_ int nMaxCount)
	{
		strx[0] = 'Z';
		char *sb = new char[4];
		strcpy_s(sb, 4, strx);
		strcpy_s((char*)lpString, nMaxCount, strx);
	}

	__declspec(dllexport) char* test6(char* strx)
	{
		strx[0] = 'Z';
		char *sb = new char[4];
		strcpy_s(sb, 4, strx);
		return sb;
	}


}