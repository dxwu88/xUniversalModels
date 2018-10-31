#include "stdafx.h"
#include "mpi.h"
#include "Constants.h"
#include "Component.h"
#include "ComponentCollection.h"
#include "Reaction.h"
#include "ReactionCollection.h"
#include "ThermoDynamics.h"
#include "ODESolver.h"

ReactorDataStruct ODESolver::m_sReactorData;

ODESolver::ODESolver()
{
	ComponentCollection* cc = new ComponentCollection();

	m_pcReactionCollection = new ReactionCollection();
}

ODESolver::~ODESolver()
{
}

void ODESolver::SetReactorData()
{
	m_sReactorData.PressureAtm = 100.0;
	m_sReactorData.ReactorTemperatureProfile = ADIABATIC;
}

bool ODESolver::Run()
{
	realtype reltol, t, tout;
	N_Vector y, abstol;
	SUNMatrix A;
	SUNLinearSolver LS;
	void *cvode_mem;
	int flag, flagr, iout;
	int rootsfound[2];

	y = abstol = NULL;
	A = NULL;
	LS = NULL;
	cvode_mem = NULL;

	/* Create serial vector of length NEQ for I.C. and abstol */
	y = N_VNew_Serial(NEQ);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
	abstol = N_VNew_Serial(NEQ);
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

	/* Initialize y */
	Ith(y, 1) = Y1;
	Ith(y, 2) = Y2;
	Ith(y, 3) = Y3;

	/* Set the scalar relative tolerance */
	reltol = RTOL;
	/* Set the vector absolute tolerance */
	Ith(abstol, 1) = ATOL1;
	Ith(abstol, 2) = ATOL2;
	Ith(abstol, 3) = ATOL3;

	/* Call CVodeCreate to create the solver memory and specify the
	* Backward Differentiation Formula and the use of a Newton iteration */
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

	/* Call CVodeInit to initialize the integrator memory and specify the
	* user's right hand side function in y'=f(t,y), the inital time T0, and
	* the initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, f, T0, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);

	// set user data
	UserDataStruct *ud = new UserDataStruct();
	ud->n = 123;
	ud->x = 234.5;

	flag = CVodeSetUserData(cvode_mem, ud);

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	* and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	/* Call CVodeRootInit to specify the root function g with 2 components */
	//flag = CVodeRootInit(cvode_mem, 2, g);
	//if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

	/* Create dense SUNMatrix for use in linear solves */
	A = SUNDenseMatrix(NEQ, NEQ);
	if (check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

	/* Create dense SUNLinearSolver object for use by CVode */
	LS = SUNDenseLinearSolver(y, A);
	if (check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	/* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	if (check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

	/* Set the user-supplied Jacobian routine Jac */
	flag = CVDlsSetJacFn(cvode_mem, Jac);
	if (check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);

	/* In loop, call CVode, print results, and test for error.
	Break out of loop when NOUT preset output times have been reached.  */
	printf(" \n3-species kinetics problem\n\n");

	iout = 0;  tout = T1;
	while (1) {
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		PrintOutput(t, Ith(y, 1), Ith(y, 2), Ith(y, 3));

		ud->n = 567;
		ud->x = 678.5;
		//flag = CVodeSetUserData(cvode_mem, ud);

		if (flag == CV_ROOT_RETURN) {
			flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
			if (check_flag(&flagr, "CVodeGetRootInfo", 1)) return(1);
			PrintRootInfo(rootsfound[0], rootsfound[1]);
		}

		if (check_flag(&flag, "CVode", 1)) break;
		if (flag == CV_SUCCESS) {
			iout++;
			tout *= TMULT;
		}

		if (iout == NOUT) break;
	}

	/* Print some final statistics */
	PrintFinalStats(cvode_mem);

	/* check the solution error */
	flag = check_ans(y, t, reltol, abstol);

	/* Free y and abstol vectors */
	N_VDestroy(y);
	N_VDestroy(abstol);

	/* Free integrator memory */
	CVodeFree(&cvode_mem);

	/* Free the linear solver memory */
	SUNLinSolFree(LS);

	/* Free the matrix memory */
	SUNMatDestroy(A);

	return(flag == 0);
}

bool ODESolver::RunA()
{
	N_Vector z;
	z = N_VNew_Serial(2);
	double yin[] = { 1.0, 2.0 };
	N_VSetArrayPointer(yin, z);
	double z0 = Ith(z, 1);
	z0 = Ith(z, 2);

	MPI_Init(NULL, NULL);

	realtype reltol, t, tout;
	N_Vector y, abstol;
	SUNMatrix A;
	SUNLinearSolver LS;
	void *cvode_mem;
	int flag, flagr, iout;
	int rootsfound[2];

	y = abstol = NULL;
	A = NULL;
	LS = NULL;
	cvode_mem = NULL;

	/* Create serial vector of length NEQ for I.C. and abstol */
	y = N_VNew_Serial(4);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
	abstol = N_VNew_Serial(4);
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

	/* Initialize y */
	Ith(y, 1) = 0.56;
	Ith(y, 2) = 1.28;
	Ith(y, 3) = 0.16;
	Ith(y, 4) = 45.0;

	/* Set the scalar relative tolerance */
	reltol = RTOL;
	/* Set the vector absolute tolerance */
	Ith(abstol, 1) = ATOL1;
	Ith(abstol, 2) = ATOL2;
	Ith(abstol, 3) = ATOL3;
	Ith(abstol, 4) = ATOL3;

	/* Call CVodeCreate to create the solver memory and specify the
	* Backward Differentiation Formula and the use of a Newton iteration */
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON); //CV_BDF
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

	/* Call CVodeInit to initialize the integrator memory and specify the
	* user's right hand side function in y'=f(t,y), the inital time T0, and
	* the initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, fa, T0, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);

	// set user data
	UserDataStruct *ud = new UserDataStruct();
	ud->n = 123;
	ud->x = 234.5;
	flag = CVodeSetUserData(cvode_mem, ud);

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	* and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	/* Call CVodeRootInit to specify the root function g with 2 components */
	//flag = CVodeRootInit(cvode_mem, 2, g);
	//if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

	/* Create dense SUNMatrix for use in linear solves */
	A = SUNDenseMatrix(4, 4);
	if (check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

	/* Create dense SUNLinearSolver object for use by CVode */
	LS = SUNDenseLinearSolver(y, A);
	if (check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	/* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	if (check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

	/* Set the user-supplied Jacobian routine Jac */
	//flag = CVDlsSetJacFn(cvode_mem, Jac);
	//if (check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);

	/* In loop, call CVode, print results, and test for error.
	Break out of loop when NOUT preset output times have been reached.  */
	printf(" \n3-species kinetics problem\n\n");

	iout = 0;  tout = 1.0;
	while (1) {
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		PrintOutputA(t, Ith(y, 1), Ith(y, 2), Ith(y, 3), Ith(y, 4));

		if (check_flag(&flag, "CVode", 1)) break;
		if (flag == CV_SUCCESS) {
			iout++;
			tout += 1.0;
		}

		if (iout == 10) break;
	}

	/* Print some final statistics */
	PrintFinalStats(cvode_mem);

	/* check the solution error */
	flag = check_ans(y, t, reltol, abstol);

	/* Free y and abstol vectors */
	N_VDestroy(y);
	N_VDestroy(abstol);

	/* Free integrator memory */
	CVodeFree(&cvode_mem);

	/* Free the linear solver memory */
	SUNLinSolFree(LS);

	/* Free the matrix memory */
	SUNMatDestroy(A);

	MPI_Finalize();

	return(flag == 0);
}

/*
*-------------------------------
* Functions called by the solver
*-------------------------------
*/

/*
* f routine. Compute function f(t,y).
*/

int ODESolver::f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype y1, y2, y3, yd1, yd3;

	UserDataStruct* ud = (UserDataStruct*)user_data;

	y1 = Ith(y, 1); 
	y2 = Ith(y, 2); 
	y3 = Ith(y, 3);

	yd1 = Ith(ydot, 1) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
	yd3 = Ith(ydot, 3) = RCONST(3.0e7)*y2*y2;
	Ith(ydot, 2) = -yd1 - yd3;

	return(0);
}

/*
* f routine. Compute function f(t,y).
*/

int ODESolver::fa(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype y1, y2, y3, y4, yd1, yd2, yd3, yd4;

	UserDataStruct* ud = (UserDataStruct*)user_data;

	double k1 = 0.950;
	//if (t > 5.0 && t <= 8.0)
	//	k1 = 0.425;
	//else if (t > 8.0 && t <= 10.0)
	//	k1 = 0.000;

	double a = 0.9082;
	double b = 1.011;
	double k2 = 3.5;
	double k3 = 18.00;
	double k4 = 37.50;
	double k5 = 1.10;

	y1 = Ith(y, 1);
	y2 = Ith(y, 2);
	y3 = Ith(y, 3);
	y4 = Ith(y, 4);

	yd1 = Ith(ydot, 1) = k1 * y1 * (1.0 - y1 / k2);
	yd2 = Ith(ydot, 2) = k3 * y1 * y4 / (k4 + y4) - a * k5 * y2;
	yd3 = Ith(ydot, 3) = k5 * y2;
	yd4 = Ith(ydot, 4) = -b * k3 * y1 * y4 / (k4 + y4);

	return(0);
}

bool ODESolver::RunReactorModel(TemperatureProfileEnum reactorTemperatureProfile, int ysize, double* yin, int tsize, double* temperatureProfile)
{
	MPI_Init(NULL, NULL);

	realtype t, tout;
	int flag, flagr, iout;
	int rootsfound[2];

	// Create serial vector of length NEQ for I.C. and abstol 
	N_Vector y;
	y = N_VNew_Serial(ysize);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return false;
	
	N_Vector abstol = NULL;
	abstol = N_VNew_Serial(ysize);
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return false;

	// Initialize y 
	N_VSetArrayPointer(yin, y);

	// Set the scalar relative tolerance
	realtype reltol = RTOL;
	
	/* Set the vector absolute tolerance */
	for (int i = 0; i < ysize; ++i)
	{
		Ith(abstol, i + 1) = ATOL3;
	}

	/* Call CVodeCreate to create the solver memory and specify the
	* Backward Differentiation Formula and the use of a Newton iteration */
	void *cvode_mem = NULL;
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON); 
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return false;

	/* Call CVodeInit to initialize the integrator memory and specify the
	* user's right hand side function in y'=f(t,y), the inital time T0, and
	* the initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, GasPhasePFRReactionDerivative, T0, y);
	if (check_flag(&flag, "CVodeInit", 1)) return false;

	// set user data
	UserDataStruct *ud = new UserDataStruct();
	ud->n = 123;
	ud->x = 234.5;
	flag = CVodeSetUserData(cvode_mem, ud);

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	* and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return false;

	/* Call CVodeRootInit to specify the root function g with 2 components */
	//flag = CVodeRootInit(cvode_mem, 2, g);
	//if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

	/* Create dense SUNMatrix for use in linear solves */
	SUNMatrix A = NULL;
	A = SUNDenseMatrix(ysize, ysize);
	if (check_flag((void *)A, "SUNDenseMatrix", 0)) return false;

	/* Create dense SUNLinearSolver object for use by CVode */
	SUNLinearSolver LS = NULL;
	LS = SUNDenseLinearSolver(y, A);
	if (check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return false;

	/* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	if (check_flag(&flag, "CVDlsSetLinearSolver", 1)) return false;

	/* Set the user-supplied Jacobian routine Jac */
	//flag = CVDlsSetJacFn(cvode_mem, Jac);
	//if (check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);

	// Controlling the number of steps the solver can take
	int maxSteps = 500;
	flag = CVodeSetMaxNumSteps(cvode_mem, maxSteps);

	/* In loop, call CVode, print results, and test for error.
	Break out of loop when NOUT preset output times have been reached.  */
	printf(" \n3-species kinetics problem\n\n");

	iout = 0;  tout = 1.0;
	while (true) {
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		PrintOutputA(t, Ith(y, 1), Ith(y, 2), Ith(y, 3), Ith(y, 4));

		if (check_flag(&flag, "CVode", 1)) break;
		if (flag == CV_SUCCESS) {
			iout++;
			tout += 1.0;
		}

		if (iout == 10) break;
	}

	/* Print some final statistics */
	PrintFinalStats(cvode_mem);

	/* check the solution error */
	flag = check_ans(y, t, reltol, abstol);

	/* Free y and abstol vectors */
	N_VDestroy(y);
	N_VDestroy(abstol);

	/* Free integrator memory */
	CVodeFree(&cvode_mem);

	/* Free the linear solver memory */
	SUNLinSolFree(LS);

	/* Free the matrix memory */
	SUNMatDestroy(A);

	MPI_Finalize();

	return (flag == 0);
}

int ODESolver::GasPhasePFRReactionDerivative(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype y1, y2, y3, y4, yd1, yd2, yd3, yd4;

	UserDataStruct* uds = (UserDataStruct*)user_data;
	ReactorDataStruct* rds = &m_sReactorData;

	//double sumOfY = y.
	// Calculate rate expressions in mols per hr per cu meter reactor
	double ConcTotal = rds->PressureAtm / (RGasM3 * rds->TempK); // R[=] m3*atm / (K*mol), ConcTotal[=] mol / m3
	double RecipRT = (1.0 / RGasKC) * (1.0 / rds->TempK - 1.0 / (rds->TRef + 273.15)); // 1 / RT(reciprocal)  expression[=](kcal / mol) - 1
	double yConc = 0.0; // ConcTotal * yMolFrac;

	double AIF = 1.0;
	typedef std::vector<Component*> ::iterator CI;
	std::vector<Component*> Components = ComponentCollection::m_vecComponents;
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		AIF += yConc * comp->m_fADS;
	}

	typedef std::vector<Reaction*> ::iterator CFR;
	typedef std::vector<std::pair<Component*, double>> ::iterator CS;

	std::vector<Reaction*> ForwardReactions = ReactionCollection::m_vecForwardReactions;
	for (CFR pr = ForwardReactions.begin(); pr != ForwardReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		int reactionIndex = reaction->m_nIndex;
		double parameter = reaction->m_fRateConstant;
		double AE = reaction->m_fActiveEnergy;

		double forwardRate = 1.0;
		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;

		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			//double gibbs = GIBBS(pair.first, temperature);
			// RXN6 : C3->C2o + C1
			//Rate(6) = FF(6) * EXP(-AE(6)*RecipRT) * (yConc(i_C3)) / (AIF)
			forwardRate *= pow(Ith(y, compIndex+1), stoich);
		}

		//Ith(ydot, reactionIndex) = parameter * exp(-AE / RecipRT) * xr / AIF;
		reaction->m_fRate = parameter * exp(-AE / RecipRT) * forwardRate / AIF;
	}

	std::vector<Reaction*> ReverseReactions = ReactionCollection::m_vecReverseReactions;
	for (CFR pr = ReverseReactions.begin(); pr != ReverseReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		int reactionIndex = reaction->m_nIndex;
		double parameter = reaction->m_fRateConstant;
		double AE = reaction->m_fActiveEnergy;

		double KEQRX = ThermoDynamics::KEQRX(reaction, rds->TempK, rds->PressureAtm);

		double reverseRate4Reactants = 1.0;
		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;
		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			// RXN2 : C3 <> C3o + H2
			// Rate(2) = FF(2) * EXP(-AE(2)*RecipRT) * (yConc(i_C3) - (yConc(i_C3o) * yConc(i_H2) / KEQRX(2, TempK, PressureAtm))) / (AIF)
			reverseRate4Reactants *= pow(Ith(y, compIndex+1), stoich);
		}

		double reverseRate4RProducts = 1.0;
		std::vector<std::pair<Component*, double>> ProductComponentStoichs = reaction->m_vecProductComponentStoich;
		for (CS pcs = ProductComponentStoichs.begin(); pcs != ProductComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;
			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			reverseRate4RProducts *= pow(Ith(y, compIndex+1), stoich);
		}

		//Ith(ydot, reactionIndex) = parameter * exp(-AE / RecipRT) * (xr - xp) / KEQRX / AIF;
		reaction->m_fRate = parameter * exp(-AE / RecipRT) * (reverseRate4Reactants - reverseRate4RProducts) / KEQRX / AIF;
	}

	// dots
	for (CFR pr = ForwardReactions.begin(); pr != ForwardReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		double rate = reaction->m_fRate;

		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;
		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			Ith(ydot, compIndex+1) -= stoich * rate;
		}
	}

	for (CFR pr = ReverseReactions.begin(); pr != ReverseReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		double rate = reaction->m_fRate;

		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;
		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			Ith(ydot, compIndex+1) -= stoich * rate;
		}

		std::vector<std::pair<Component*, double>> ProductComponentStoichs = reaction->m_vecProductComponentStoich;
		for (CS pcs = ProductComponentStoichs.begin(); pcs != ProductComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;
			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			Ith(ydot, compIndex+1) += stoich * rate;
		}
	}

	for (int i = 0; i < Components.size(); ++i)
	{
		Ith(ydot, i + 1) /= ConcTotal;
	}

	if (rds->ReactorTemperatureProfile == ADIABATIC)
	{
		// use ideal gas enthalpy to calculate the adiabatic temperature change, KJ / hr
		double DeltaEnthalpy = ThermoDynamics::IdealGasEnthalpy((double*)ydot, rds->TempK);

		// DeltaH put on molar basis, KJ/Kg-mol
		double TotalMoles = 1.0;
		double DeltaH = DeltaEnthalpy / TotalMoles;
		// Cp is in KJ / (Kg - mol * K)
		double Cp = ThermoDynamics::HeatCapacity((double*)y, rds->TempK);

		Ith(ydot, Components.size() + 1) = -DeltaH / Cp;
	}

	return 0;
}

int ODESolver::GasAndLiquidPhasePFRReactionDerivative(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	UserDataStruct* uds = (UserDataStruct*)user_data;
	ReactorDataStruct* rds = &m_sReactorData;

	double Q_L = 0;
	double VaporFraction = 0.0;
	//rxns.CalculateFlash(y, sum, VaporFraction, SumMolesG, SumMolesL, Q_L); // Calculating Flash and separating streams																			   //rxns.VFvector.push_back(VaporFraction);
	//m_pcReactionCollection->SetVaporFraction(VaporFraction);

	// Calculate rate expressions in mols per hr per cu meter reactor
	double ConcTotal = rds->PressureAtm / (RGasM3 * rds->TempK); // R[=] m3*atm / (K*mol), ConcTotal[=] mol / m3
	double RecipRT = (1.0 / RGasKC) * (1.0 / rds->TempK - 1.0 / (rds->TRef + 273.15)); // 1 / RT(reciprocal)  expression[=](kcal / mol) - 1
	double yConc = 0.0; // ConcTotal * yMolFrac;

	double AIF = 1.0;
	typedef std::vector<Component*> ::iterator CI;
	std::vector<Component*> Components = ComponentCollection::m_vecComponents;
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		AIF += yConc * comp->m_fADS;
	}

	typedef std::vector<Reaction*> ::iterator CFR;
	typedef std::vector<std::pair<Component*, double>> ::iterator CS;

	std::vector<Reaction*> ForwardReactions = ReactionCollection::m_vecForwardReactions;
	for (CFR pr = ForwardReactions.begin(); pr != ForwardReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		int reactionIndex = reaction->m_nIndex;
		double parameter = reaction->m_fRateConstant;
		double AE = reaction->m_fActiveEnergy;

		double forwardRate = 1.0;
		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;

		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			//double gibbs = GIBBS(pair.first, temperature);
			// RXN6 : C3->C2o + C1
			//Rate(6) = FF(6) * EXP(-AE(6)*RecipRT) * (yConc(i_C3)) / (AIF)
			forwardRate *= pow(Ith(y, compIndex + 1), stoich);
		}

		//Ith(ydot, reactionIndex) = parameter * exp(-AE / RecipRT) * xr / AIF;
		reaction->m_fRate = parameter * exp(-AE / RecipRT) * forwardRate / AIF;
	}

	std::vector<Reaction*> ReverseReactions = ReactionCollection::m_vecReverseReactions;
	for (CFR pr = ReverseReactions.begin(); pr != ReverseReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		int reactionIndex = reaction->m_nIndex;
		double parameter = reaction->m_fRateConstant;
		double AE = reaction->m_fActiveEnergy;

		double KEQRX = ThermoDynamics::KEQRX(reaction, rds->TempK, rds->PressureAtm);

		double reverseRate4Reactants = 1.0;
		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;
		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			// RXN2 : C3 <> C3o + H2
			// Rate(2) = FF(2) * EXP(-AE(2)*RecipRT) * (yConc(i_C3) - (yConc(i_C3o) * yConc(i_H2) / KEQRX(2, TempK, PressureAtm))) / (AIF)
			reverseRate4Reactants *= pow(Ith(y, compIndex + 1), stoich);
		}

		double reverseRate4RProducts = 1.0;
		std::vector<std::pair<Component*, double>> ProductComponentStoichs = reaction->m_vecProductComponentStoich;
		for (CS pcs = ProductComponentStoichs.begin(); pcs != ProductComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;
			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			reverseRate4RProducts *= pow(Ith(y, compIndex + 1), stoich);
		}

		//Ith(ydot, reactionIndex) = parameter * exp(-AE / RecipRT) * (xr - xp) / KEQRX / AIF;
		reaction->m_fRate = parameter * exp(-AE / RecipRT) * (reverseRate4Reactants - reverseRate4RProducts) / KEQRX / AIF;
	}

	// dots
	for (CFR pr = ForwardReactions.begin(); pr != ForwardReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		double rate = reaction->m_fRate;

		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;
		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			Ith(ydot, compIndex + 1) -= stoich * rate;
		}
	}

	for (CFR pr = ReverseReactions.begin(); pr != ReverseReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		double rate = reaction->m_fRate;

		std::vector<std::pair<Component*, double>> ReactantComponentStoichs = reaction->m_vecReactantComponentStoich;
		for (CS pcs = ReactantComponentStoichs.begin(); pcs != ReactantComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			Ith(ydot, compIndex + 1) -= stoich * rate;
		}

		std::vector<std::pair<Component*, double>> ProductComponentStoichs = reaction->m_vecProductComponentStoich;
		for (CS pcs = ProductComponentStoichs.begin(); pcs != ProductComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;
			Component* comp = pair.first;
			int compIndex = comp->m_nIndex;
			double stoich = pair.second;

			Ith(ydot, compIndex + 1) += stoich * rate;
		}
	}

	for (int i = 0; i < Components.size(); ++i)
	{
		Ith(ydot, i + 1) /= ConcTotal;
	}

	if (rds->ReactorTemperatureProfile == ADIABATIC)
	{
		// use ideal gas enthalpy to calculate the adiabatic temperature change, KJ / hr
		double DeltaEnthalpy = ThermoDynamics::IdealGasEnthalpy((double*)ydot, rds->TempK);

		// DeltaH put on molar basis, KJ/Kg-mol
		double TotalMoles = 1.0;
		double DeltaH = DeltaEnthalpy / TotalMoles;
		// Cp is in KJ / (Kg - mol * K)
		double Cp = ThermoDynamics::HeatCapacity((double*)y, rds->TempK);

		Ith(ydot, Components.size() + 1) = -DeltaH / Cp;
	}

	return 0;
}

/***
	subroutine CalcRates()
		!*****************************************************************
		!This subroutine calculates the reaction rates at the current
		!conditions for the integration.
		use ReactionData, only: PressureAtm, &
		TempK, &
		TRef, &
		yMolFrac, &
		FF, &
		AE, &
		Rate, &
		TotalComps, &
		TotalRxns, &
		ADS

		use CompData
		use ThermoRoutines
		implicit none

		integer i

		real * 8 RGasKC
		real * 8 RecipRT
		real * 8 ConcTotal
		real * 8 AIF
		real * 8 yConc(TotalComps)

		!Calculate rate expressions in mols per hr per cu meter reactor
		ConcTotal = PressureAtm / (RGasM3 * TempK)                                 !R[=] m3*atm / (K*mol), ConcTotal[=] mol / m3
		RGasKC = 1.9872041D - 03                                                      !R[=] kcal / (K*mol)
		RecipRT = (1.0D0 / RGasKC) * (1.0D0 / TempK - 1.0D0 / (TRef + 273.15D0))          !1 / RT(reciprocal)  expression[=](kcal / mol) - 1
		yConc = ConcTotal * yMolFrac

		AIF = 1.0d0
		do i = 1, TotalComps
			AIF = AIF + yConc(i)*ADS(i)
			enddo
			AIF = AIF * *2

			!Reaction rate equations
			Rate = 0.0D0

			!RXN1 :  C2 <> C2o + H2
			Rate(1) = FF(1) * EXP(-AE(1)*RecipRT) * (yConc(i_C2) - (yConc(i_C2o) * yConc(i_H2) / KEQRX(1, TempK, PressureAtm))) / (AIF)

			!RXN2 : C3 <> C3o + H2
			Rate(2) = FF(2) * EXP(-AE(2)*RecipRT) * (yConc(i_C3) - (yConc(i_C3o) * yConc(i_H2) / KEQRX(2, TempK, PressureAtm))) / (AIF)

			!RXN3 : nC4 <> nC4o + H2
			Rate(3) = FF(3) * EXP(-AE(3)*RecipRT) * (yConc(i_nC4) - (yConc(i_nC4o) * yConc(i_H2) / KEQRX(3, TempK, PressureAtm))) / (AIF)

			!RXN4 : nC5 <> nC5o + H2
			Rate(4) = FF(4) * EXP(-AE(4)*RecipRT) * (yConc(i_nC5) - (yConc(i_nC5o) * yConc(i_H2) / KEQRX(4, TempK, PressureAtm))) / (AIF)

			!RXN5 : nC6 <> nC6o + H2
			Rate(5) = FF(5) * EXP(-AE(5)*RecipRT) * (yConc(i_nC6) - (yConc(i_nC6o) * yConc(i_H2) / KEQRX(5, TempK, PressureAtm))) / (AIF)

			!RXN6 : C3->C2o + C1
			Rate(6) = FF(6) * EXP(-AE(6)*RecipRT) * (yConc(i_C3)) / (AIF)

			!RXN7 : nC4->C3o + C1
			Rate(7) = FF(7) * EXP(-AE(7)*RecipRT) * (yConc(i_nC4)) / (AIF)

			!RXN8 : nC4->C2o + C2
			Rate(8) = FF(8) * EXP(-AE(8)*RecipRT) * (yConc(i_nC4)) / (AIF)

			!RXN9 : nC5->nC4o + C1
			Rate(9) = FF(9) * EXP(-AE(9)*RecipRT) * (yConc(i_nC5)) / (AIF)

			!RXN10 : nC5->C3o + C2
			Rate(10) = FF(10) * EXP(-AE(10)*RecipRT) * (yConc(i_nC5)) / (AIF)

			!RXN11 : nC5->C2o + C3
			Rate(11) = FF(11) * EXP(-AE(11)*RecipRT) * (yConc(i_nC5)) / (AIF)

			!RXN12 : nC6o <> MCC5
			Rate(12) = FF(12) * EXP(-AE(12)*RecipRT) * (yConc(i_nC6o) - (yConc(i_MCC5) / KEQRX(12, TempK, PressureAtm))) / (AIF)

			!RXN13 : C3o + C3 <> iC6
			Rate(13) = FF(13) * EXP(-AE(13)*RecipRT) * (yConc(i_C3o) * yConc(i_C3) - (yConc(i_iC6) / KEQRX(13, TempK, PressureAtm))) / (AIF)

			!RXN14 : nC6 <> iC6
			Rate(14) = FF(14) * EXP(-AE(14)*RecipRT) * (yConc(i_nC6) - (yConc(i_iC6) / KEQRX(14, TempK, PressureAtm))) / (AIF)

			!RXN15 : C3o + nC4 <> iC7
			Rate(15) = FF(15) * EXP(-AE(15)*RecipRT) * (yConc(i_C3o) * yConc(i_nC4) - (yConc(i_iC7) / KEQRX(15, TempK, PressureAtm))) / (AIF)

			!RXN16 : 2 * C3o <> nC6o
			Rate(16) = FF(16) * EXP(-AE(16)*RecipRT) * (yConc(i_C3o)**2 - (yConc(i_nC6o) / KEQRX(16, TempK, PressureAtm))) / (AIF)

			!RXN17 : 2 * nC4o <> nC5o + C3o
			Rate(17) = FF(17) * EXP(-AE(17)*RecipRT) * (yConc(i_nC4o)**2 - (yConc(i_nC5o) * yConc(i_C3o) / KEQRX(17, TempK, PressureAtm))) / (AIF)

			!RXN18 : 2 * nC4o <> nC6o + C2o
			Rate(18) = FF(18) * EXP(-AE(18)*RecipRT) * (yConc(i_nC4o)**2 - (yConc(i_nC6o) * yConc(i_C2o) / KEQRX(18, TempK, PressureAtm))) / (AIF)

			!RXN19 : nC4o + C3o <> nC5o + C2o
			Rate(19) = FF(19) * EXP(-AE(19)*RecipRT) * (yConc(i_nC4o) * yConc(i_C3o) - (yConc(i_nC5o) * yConc(i_C2o) / KEQRX(19, TempK, PressureAtm))) / (AIF)

			!RXN20 : nC4o <> 2 * C2o
			Rate(20) = FF(20) * EXP(-AE(20)*RecipRT) * (yConc(i_nC4o) - (yConc(i_C2o)**2 / KEQRX(20, TempK, PressureAtm))) / (AIF)

			!RXN21 : nC5o <> C3o + C2o
			Rate(21) = FF(21) * EXP(-AE(21)*RecipRT) * (yConc(i_nC5o) - (yConc(i_C3o) * yConc(i_C2o) / KEQRX(21, TempK, PressureAtm))) / (AIF)

			!RXN22 : nC6o <> nC4o + C2o
			Rate(22) = FF(22) * EXP(-AE(22)*RecipRT) * (yConc(i_nC6o) - (yConc(i_nC4o) * yConc(i_C2o) / KEQRX(22, TempK, PressureAtm))) / (AIF)

			!RXN23 : nC6o + C2o <> nC5o + C3o
			Rate(23) = FF(23) * EXP(-AE(23)*RecipRT) * (yConc(i_nC6o) * yConc(i_C2o) - (yConc(i_nC5o) * yConc(i_C3o) / KEQRX(23, TempK, PressureAtm))) / (AIF)

			!RXN24 : nC6o <> A6 + 3 * H2
			Rate(24) = FF(24) * EXP(-AE(24)*RecipRT) * (yConc(i_nC6o) - (yConc(i_A6) * yConc(i_H2)**3 / KEQRX(24, TempK, PressureAtm))) / (AIF)

			!RXN25 : nC4o + C2o <> A6 + 3 * H2
			Rate(25) = FF(25) * EXP(-AE(25)*RecipRT) * (yConc(i_nC4o) * yConc(i_C2o) - (yConc(i_A6) * yConc(i_H2)**3 / KEQRX(25, TempK, PressureAtm))) / (AIF)

			!RXN26 : C3o + nC4o <> A7 + 3 * H2
			Rate(26) = FF(26) * EXP(-AE(26)*RecipRT) * (yConc(i_C3o) * yConc(i_nC4o) - (yConc(i_A7) * yConc(i_H2)**3 / KEQRX(26, TempK, PressureAtm))) / (AIF)

			!RXN27 : C2o + nC5o <> A7 + 3 * H2
			Rate(27) = FF(27) * EXP(-AE(27)*RecipRT) * (yConc(i_C2o) * yConc(i_nC5o) - (yConc(i_A7) * yConc(i_H2)**3 / KEQRX(27, TempK, PressureAtm))) / (AIF)

			!RXN28 : C2o + nC6o <> A8 + 3 * H2
			Rate(28) = FF(28) * EXP(-AE(28)*RecipRT) * (yConc(i_C2o) * yConc(i_nC6o) - (yConc(i_A8) * yConc(i_H2)**3 / KEQRX(28, TempK, PressureAtm))) / (AIF)

			!RXN29 : C3o + nC5o <> A8 + 3 * H2
			Rate(29) = FF(29) * EXP(-AE(29)*RecipRT) * (yConc(i_C3o) * yConc(i_nC5o) - (yConc(i_A8) * yConc(i_H2)**3 / KEQRX(29, TempK, PressureAtm))) / (AIF)

			!RXN30 : 2 * nC4o <> A8 + 3 * H2
			Rate(30) = FF(30) * EXP(-AE(30)*RecipRT) * (yConc(i_nC4o)**2 - (yConc(i_A8) * yConc(i_H2)**3 / KEQRX(30, TempK, PressureAtm))) / (AIF)

			!RXN31 : 2 * A10 <> A11 + A9
			Rate(31) = FF(31) * EXP(-AE(31)*RecipRT) * (yConc(i_A10)**2 - (yConc(i_A11) * yConc(i_A9) / KEQRX(31, TempK, PressureAtm))) / (AIF)

			!RXN32 : A7 + A8 <> A6 + A9
			Rate(32) = FF(32) * EXP(-AE(32)*RecipRT) * (yConc(i_A7) * yConc(i_A8) - (yConc(i_A6) * yConc(i_A9) / KEQRX(32, TempK, PressureAtm))) / (AIF)

			!RXN33 : A8 + A9 <> A7 + A10
			Rate(33) = FF(33) * EXP(-AE(33)*RecipRT) * (yConc(i_A8) * yConc(i_A9) - (yConc(i_A7) * yConc(i_A10) / KEQRX(33, TempK, PressureAtm))) / (AIF)

			!RXN34 : A8 + A10 <> A7 + A11
			Rate(34) = FF(34) * EXP(-AE(34)*RecipRT) * (yConc(i_A8) * yConc(i_A10) - (yConc(i_A7) * yConc(i_A11) / KEQRX(34, TempK, PressureAtm))) / (AIF)

			!RXN35 : A8 + A11 <> A7 + A12
			Rate(35) = FF(35) * EXP(-AE(35)*RecipRT) * (yConc(i_A8) * yConc(i_A11) - (yConc(i_A7) * yConc(i_A12) / KEQRX(35, TempK, PressureAtm))) / (AIF)

			!RXN36 : C3o + A6 <> INDANE + H2
			Rate(36) = FF(36) * EXP(-AE(36)*RecipRT) * (yConc(i_C3o) * yConc(i_A6) - (yConc(i_INDANE) * yConc(i_H2) / KEQRX(36, TempK, PressureAtm))) / (AIF)

			!RXN37 : C3o + A6 <> INDENE + 2 * H2
			Rate(37) = FF(37) * EXP(-AE(37)*RecipRT) * (yConc(i_C3o) * yConc(i_A6) - (yConc(i_INDENE) * yConc(i_H2)**2 / KEQRX(37, TempK, PressureAtm))) / (AIF)

			!RXN38 : nC4o + A6 <> Naph + 3 * H2
			Rate(38) = FF(38) * EXP(-AE(38)*RecipRT) * (yConc(i_nC4o) * yConc(i_A6) - (yConc(i_Naph) * yConc(i_H2)**3 / KEQRX(38, TempK, PressureAtm))) / (AIF)

			!RXN39 : nC4o + A7 <> MNaph + 3 * H2
			Rate(39) = FF(39) * EXP(-AE(39)*RecipRT) * (yConc(i_nC4o) * yConc(i_A7) - (yConc(i_MNaph) * yConc(i_H2)**3 / KEQRX(39, TempK, PressureAtm))) / (AIF)

			!RXN40 : nC5o + A7 <> MMNaph + 3 * H2
			Rate(40) = FF(40) * EXP(-AE(40)*RecipRT) * (yConc(i_nC5o) * yConc(i_A7) - (yConc(i_MMNaph) * yConc(i_H2)**3 / KEQRX(40, TempK, PressureAtm))) / (AIF)

			!RXN41 : C3o + Naph <> Fluorene + 2 * H2
			Rate(41) = FF(41) * EXP(-AE(41)*RecipRT) * (yConc(i_C3o) * yConc(i_Naph) - (yConc(i_Fluorene) * yConc(i_H2)**2 / KEQRX(41, TempK, PressureAtm))) / (AIF)

			!RXN42 : C3o + MNaph <> Phen + 3 * H2
			Rate(42) = FF(42) * EXP(-AE(42)*RecipRT) * (yConc(i_C3o) * yConc(i_MNaph) - (yConc(i_Phen) * yConc(i_H2)**3 / KEQRX(42, TempK, PressureAtm))) / (AIF)

			!RXN43 : nC4o + MNaph <> MAnth + 3 * H2
			Rate(43) = FF(43) * EXP(-AE(43)*RecipRT) * (yConc(i_nC4o) * yConc(i_MNaph) - (yConc(i_MAnth) * yConc(i_H2)**3 / KEQRX(43, TempK, PressureAtm))) / (AIF)


			end subroutine CalcRates
***/

/*
* g routine. Compute functions g_i(t,y) for i = 0,1.
*/

int ODESolver::g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
	realtype y1, y3;

	y1 = Ith(y, 1); y3 = Ith(y, 3);
	gout[0] = y1 - RCONST(0.0001);
	gout[1] = y3 - RCONST(0.01);

	return(0);
}

/*
* Jacobian routine. Compute J(t,y) = df/dy. *
*/

int ODESolver::Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
	void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	realtype y1, y2, y3;

	y1 = Ith(y, 1); y2 = Ith(y, 2); y3 = Ith(y, 3);

	IJth(J, 1, 1) = RCONST(-0.04);
	IJth(J, 1, 2) = RCONST(1.0e4)*y3;
	IJth(J, 1, 3) = RCONST(1.0e4)*y2;

	IJth(J, 2, 1) = RCONST(0.04);
	IJth(J, 2, 2) = RCONST(-1.0e4)*y3 - RCONST(6.0e7)*y2;
	IJth(J, 2, 3) = RCONST(-1.0e4)*y2;

	IJth(J, 3, 1) = ZERO;
	IJth(J, 3, 2) = RCONST(6.0e7)*y2;
	IJth(J, 3, 3) = ZERO;

	return(0);
}

/*
*-------------------------------
* Private helper functions
*-------------------------------
*/

void ODESolver::PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
	printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
	printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
	printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

	return;
}

void ODESolver::PrintOutputA(realtype t, realtype y1, realtype y2, realtype y3, realtype y4)
{
	printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3, y4);

	return;
}

void ODESolver::PrintRootInfo(int root_f1, int root_f2)
{
	printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);

	return;
}

/*
* Get and print some final statistics
*/

void ODESolver::PrintFinalStats(void *cvode_mem)
{
	long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
	int flag;

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	check_flag(&flag, "CVodeGetNumRhsEvals", 1);
	flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
	check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	check_flag(&flag, "CVodeGetNumErrTestFails", 1);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

	flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
	check_flag(&flag, "CVDlsGetNumJacEvals", 1);
	flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
	check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

	flag = CVodeGetNumGEvals(cvode_mem, &nge);
	check_flag(&flag, "CVodeGetNumGEvals", 1);

	printf("\nFinal Statistics:\n");
	printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
		nst, nfe, nsetups, nfeLS, nje);
	printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
		nni, ncfn, netf, nge);
}

/*
* Check function return value...
*   opt == 0 means SUNDIALS function allocates memory so check if
*            returned NULL pointer
*   opt == 1 means SUNDIALS function returns a flag so check if
*            flag >= 0
*   opt == 2 means function allocates memory so check if returned
*            NULL pointer
*/

int ODESolver::check_flag(void *flagvalue, const char *funcname, int opt)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
			funcname);
		return(1);
	}

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *)flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
				funcname, *errflag);
			return(1);
		}
	}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
			funcname);
		return(1);
	}

	return(0);
}

/* compare the solution at the final time 4e10s to a reference solution computed
using a relative tolerance of 1e-8 and absoltue tolerance of 1e-14 */
int ODESolver::check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol)
{
	int      passfail = 0;        /* answer pass (0) or fail (1) flag */
	N_Vector ref;               /* reference solution vector        */
	N_Vector ewt;               /* error weight vector              */
	realtype err;               /* wrms error                       */
	realtype ONE = RCONST(1.0);

	/* create reference solution and error weight vectors */
	ref = N_VClone(y);
	ewt = N_VClone(y);

	/* set the reference solution data */
	NV_Ith_S(ref, 0) = RCONST(5.2083495894337328e-08);
	NV_Ith_S(ref, 1) = RCONST(2.0833399429795671e-13);
	NV_Ith_S(ref, 2) = RCONST(9.9999994791629776e-01);

	/* compute the error weight vector, loosen atol */
	N_VAbs(ref, ewt);
	N_VLinearSum(rtol, ewt, RCONST(10.0), atol, ewt);
	if (N_VMin(ewt) <= ZERO) {
		fprintf(stderr, "\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n\n");
		return(-1);
	}
	N_VInv(ewt, ewt);

	/* compute the solution error */
	N_VLinearSum(ONE, y, -ONE, ref, ref);
	err = N_VWrmsNorm(ref, ewt);

	/* is the solution within the tolerances? */
	passfail = (err < ONE) ? 0 : 1;

	if (passfail) {
		fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%g \n\n", err);
	}

	/* Free vectors */
	N_VDestroy(ref);
	N_VDestroy(ewt);

	return(passfail);
}
