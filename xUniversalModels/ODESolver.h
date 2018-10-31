#pragma once

/* -----------------------------------------------------------------
* Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
*                Radu Serban @ LLNL
* -----------------------------------------------------------------
* Example problem :
*
* The following is a simple example problem, with the coding
* needed for its solution by CVODE.The problem is from
* chemical kinetics, and consists of the following three rate
* equations :
	*	dy1 / dt = -.04*y1 + 1.e4*y2*y3
	*	dy2 / dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2) ^ 2
	*	dy3 / dt = 3.e7*(y2) ^ 2
	* on the interval from t = 0.0 to t = 4.e10, with initial
	* conditions : y1 = 1.0, y2 = y3 = 0. The problem is stiff.
	* While integrating the system, we also use the rootfinding
	* feature to find the points at which y1 = 1e-4 or at which
	* y3 = 0.01.This program solves the problem with the BDF method,
	* Newton iteration with the SUNDENSE dense linear solver, and a
	* user - supplied Jacobian routine.
	* It uses a scalar relative tolerance and a vector absolute
	* tolerance.Output is printed in decades from t = .4 to t = 4.e10.
	* Run statistics(optional outputs) are printed at the end.
	* ---------------------------------------------------------------- - */

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "ReactionCollection.h"

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */

#define NEQ   3                /* number of equations  */
#define Y1    RCONST(1.0)      /* initial y components */
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-8)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-14)
#define ATOL3 RCONST(1.0e-6)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.4)      /* first output time      */
#define TMULT RCONST(10.0)     /* output time factor     */
#define NOUT  12               /* number of output times */

#include "StructAndEnum.h"

class ODESolver
{
public:
	static ReactorDataStruct m_sReactorData;

private:
	ReactionCollection* m_pcReactionCollection;

public:
	ODESolver();
	~ODESolver();

	void SetReactorData();

	bool Run();
	bool RunA();
	bool RunReactorModel(TemperatureProfileEnum reactorTemperatureProfile, int ysize, double* yin, int tsize, double* temperatureProfile);

	/* Functions Called by the Solver */

	static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
	static int fa(realtype t, N_Vector y, N_Vector ydot, void *user_data);
	static int GasPhasePFRReactionDerivative(realtype t, N_Vector y, N_Vector ydot, void *user_data);
	static int GasAndLiquidPhasePFRReactionDerivative(realtype t, N_Vector y, N_Vector ydot, void *user_data);

	static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

	static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
	void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

	/* Private functions to output results */

	static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);
	static void PrintOutputA(realtype t, realtype y1, realtype y2, realtype y3, realtype y4);
	static void PrintRootInfo(int root_f1, int root_f2);

	/* Private function to print final statistics */

	static void PrintFinalStats(void *cvode_mem);

	/* Private function to check function return values */

	static int check_flag(void *flagvalue, const char *funcname, int opt);

	/* Private function to check computed solution */

	static int check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol);


};

