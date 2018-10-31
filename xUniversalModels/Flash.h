#pragma once
class Flash
{
public:
	Flash();
	~Flash();

private:
	// VLE Section
	double *speciesCP_l;
	double *speciesH_l;
	double *speciesG_l;
	double *deltaGrxn_l;
	double *deltaHrxn_l;
	double *Pvap;
	double *Moles_L;
	double *Moles_G;
	double *MWi;
	double *Densityi;

public:
	static const double DELTAS[];

private:
	void CalcPengRobinsonPvaps(vector<double> &Pvaps);
	double LeeKeslerPvap(double T, double Tc, double Pc, double w);
	double NewtonRaphsonBisectionForRachfordRice(double pressure, N_Vector y, double sum, double xacc, int maxiter);

public:
	bool CalculateFlash(N_Vector y, double temperature, double pressure, double SumMoles, double &VaporFraction, double &SumMolesG, double &SumMolesL, double &Q_L);

};

