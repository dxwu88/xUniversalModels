// TestSundials.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
//#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include "mpi.h"

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

int main(int argc, char *argv[])
{
	ODESolver *os = new ODESolver();
	os->RunA();

	int N = 200;
	realtype t0 = 0;
	realtype Tfinal = 10;
	realtype theta0 = 0.0; // atof(argv[1]);
	realtype reltol = 1e-6;
	realtype abstol = 1e-8;
	realtype t;
	int flag, k;
	N_Vector y = NULL;
	void* cvode_mem = NULL;
	// Create serial vector of length NEQ for I.C. 
	y = N_VNew_Serial(2);

	NV_Ith_S(y, 0) = 0.1; // theta0;
	NV_Ith_S(y, 1) = 0;
	
	N_Vector abstolV = N_VNew_Serial(2);
	NV_Ith_S(abstolV, 0) = 1e-8;
	NV_Ith_S(abstolV, 1) = 1e-8;

	//MPI::Init(&argc, &argv);
	MPI_Init(NULL, NULL);

	// Set up solver 
	cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
	if (cvode_mem == 0) {
		fprintf(stderr, "Error in CVodeMalloc: could not allocate\n");
		return -1;
	}
	// Call CVodeMalloc to initialize the integrator memory 
	flag = CVodeInit(cvode_mem, f, t0, y);
	if (flag < 0) {
		fprintf(stderr, "Error in CVodeMalloc: %d\n", flag);
		return -1;
	}

	flag = CVodeSVtolerances(cvode_mem, reltol, abstolV);

	// In loop, call CVode, print results, and test for error. 
	for (k = 1; k < N; ++k) {
		realtype tout = k * Tfinal / N;
		if (CVode(cvode_mem, tout, y, &t, CV_NORMAL) < 0) {
			fprintf(stderr, "Error in CVode: %d\n", flag);
			return -1;
		}
		printf("%g %.16e %.16e\n", t, NV_Ith_S(y, 0), NV_Ith_S(y, 1));
	}
	N_VDestroy_Serial(y); // Free y vector 
	CVodeFree(&cvode_mem); // Free integrator memory 
    
	MPI_Finalize();

	return 0;
}
