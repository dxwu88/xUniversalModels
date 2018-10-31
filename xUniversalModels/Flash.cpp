#include "stdafx.h"
#include "Component.h"
#include "ComponentCollection.h"
#include "ThermoDynamics.h"
#include "Flash.h"

const double Flash::DELTAS[] = { 0.5, 1.0, 4.0, 7.0, 10.0, 10.0, 20.0, 20.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0 };

Flash::Flash()
{
}

Flash::~Flash()
{
}

bool Flash::CalculateFlash(N_Vector y, double temperature, double pressure, double SumMoles, double &VaporFraction, double &SumMolesG, double &SumMolesL, double &Q_L)
{
	std::vector<Component*> Components = ComponentCollection::m_vecComponents;

	// Calculating for the fraction of the feed that is vaporized
	VaporFraction = NewtonRaphsonBisectionForRachfordRice(pressure, y, SumMoles, 0.00001, 100); 

	if (VaporFraction >= 0.0 && VaporFraction <= 1.0)
	{
		double MoleFractionsG = 0.0;
		double MoleFractionsL = 0.0;

		// Different reactors matter here, need to change this to account for different reactors this is coming from
		// For PFR mols/s differentiating by volume
		// For Batch mols/L differentiating by time
		/***
		for (int i = 0; i < NumMolecules; i++)
		{
			MoleFractionsL = (Ith(y, i + 1) / SumMoles) / (1 + ((Pvap[i] / pressure)*VaporFraction));
			MoleFractionsG = (Pvap[i] / pressure) * MoleFractionsL;
			//Use Henry's Law to calculate H2
			//Data from http://www.henrys-law.org/henry-3.0.pdf
			if (Components[i]->m_strName.compare("Hydrogen") == 0)
			{
				Moles_L[i] = 0.00078*exp(500.0 * ((1.0 / temperature) - (1.0 / 298.15))) * pressure / SumMoles; //Pres in atm. k_H units are M/atm
				SumMolesL += Moles_L[i];
				Q_L += (Moles_L[i] * MWi[i]) / (Densityi[i] * 1000.0);
				Moles_G[i] = Ith(y, i + 1) - Moles_L[i];
				SumMolesG += Moles_G[i];
			}
			else
			{
				Moles_L[i] = MoleFractionsL * SumMoles * (1.0 - VaporFraction); // Returning the Liquid Phase
				SumMolesL += MoleFractionsL * SumMoles * (1.0 - VaporFraction); // Adding up the moles in the liquid phase
				Q_L += (Moles_L[i] * MWi[i]) / (Densityi[i] * 1000.0);
				Moles_G[i] = MoleFractionsG * SumMoles * VaporFraction;
				SumMolesG += MoleFractionsG * SumMoles * VaporFraction;
			}
		}
		***/

		typedef std::vector<Component*> ::iterator CI;
		for (CI p = Components.begin(); p != Components.end(); ++p)
		{
			Component* comp = *p;
			
			int i = comp->m_nIndex;

			MoleFractionsL = (Ith(y, i + 1) / SumMoles) / (1 + ((Pvap[i] / pressure)*VaporFraction));
			MoleFractionsG = (Pvap[i] / pressure) * MoleFractionsL;

			if (comp->m_strName.compare("Hydrogen") == 0)
			{
				Moles_L[i] = 0.00078*exp(500.0 * ((1.0 / temperature) - (1.0 / 298.15))) * pressure / SumMoles; //Pres in atm. k_H units are M/atm
				SumMolesL += Moles_L[i];
				Q_L += (Moles_L[i] * MWi[i]) / (Densityi[i] * 1000.0);
				Moles_G[i] = Ith(y, i + 1) - Moles_L[i];
				SumMolesG += Moles_G[i];
			}
			else
			{
				Moles_L[i] = MoleFractionsL * SumMoles * (1.0 - VaporFraction); // Returning the Liquid Phase
				SumMolesL += MoleFractionsL * SumMoles * (1.0 - VaporFraction); // Adding up the moles in the liquid phase
				Q_L += (Moles_L[i] * MWi[i]) / (Densityi[i] * 1000.0);
				Moles_G[i] = MoleFractionsG * SumMoles * VaporFraction;
				SumMolesG += MoleFractionsG * SumMoles * VaporFraction;
			}
		}
		
		return true;
	}

	return false;
}

//This function uses the Newton-Raphson method with bisection to solve the Rachford-Rice Equation 
//P is the system pressure, xacc is the accuracy required, maxiter is the maximum number of iterations of Newton-Raphson method
double Flash::NewtonRaphsonBisectionForRachfordRice(double pressure, N_Vector y, double sum, double xacc, int maxiter)
{
	typedef std::vector<Component*> ::iterator CI;

	std::vector<Component*> Components = ComponentCollection::m_vecComponents;

	int i, j; //Loop iterators
	double dx, dxold, fh = 0.0, fl = 0.0, f = 0.0, df = 0.0; //f is the function value and df is the value of the derivative of f
	double temp, xh, xl;
	double vf; //The vapor fraction (answer))
	double ub = 1.0; //Upper bound of the solution
	double lb = 0.0; //Lower bound of the solution

	//Evaluate the Rachford-Rice Equation at the minimum and maximum points. Pvap does not exist (=-1) for reaction families, which are species as defined currently
	//for (i = 0; i < NumMolecules; i++) {
	//	fl += (Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1) / (1 + lb * ((Pvap[i] / pressure) - 1));
	//	fh += (Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1) / (1 + ub * ((Pvap[i] / pressure) - 1));
	//}
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		fl += (Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1) / (1 + lb * ((Pvap[i] / pressure) - 1));
		fh += (Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1) / (1 + ub * ((Pvap[i] / pressure) - 1));
	}

	//Return the lower or upper bound if the function value is within the error xacc at those points
	if (fabs(fl) < xacc) return lb;
	if (fabs(fh) < xacc) return ub;
	//Orient the search so that f(xl) < 0.
	if (fl < 0.0) {
		xl = lb;
		xh = ub;
	}
	else {
		xh = lb;
		xl = ub;
	}
	vf = 0.9; //Initialize the guess for vapor fraction,
	dxold = fabs(ub - lb); //the “stepsize before last”
	dx = dxold; //and the last step.

				//Calculate Rachford-Rice Equation value at guess. Pvap does not exist (=-1) for reaction families, which are species as defined currently
	//for (i = 0; i< NumMolecules; i++) {
	//	f += (Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1) / (1 + vf * ((Pvap[i] / pressure) - 1));
	//	df += -(Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1)*((Pvap[i] / pressure) - 1) / ((1 + vf * ((Pvap[i] / pressure) - 1))*(1 + vf * ((Pvap[i] / pressure) - 1)));
	//}
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		f += (Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1) / (1 + vf * ((Pvap[i] / pressure) - 1));
		df += -(Ith(y, i + 1) / sum) * ((Pvap[i] / pressure) - 1)*((Pvap[i] / pressure) - 1) / ((1 + vf * ((Pvap[i] / pressure) - 1))*(1 + vf * ((Pvap[i] / pressure) - 1)));
	}

	//Newton-Raphson method with bisection from NUMERICAL RECIPES IN C
	//Loop over allowed iterations.
	for (i = 0; i < maxiter; i++) {
		//Bisect if Newton out of range or not decreasing fast enough.
		if ((((vf - xh)*df - f)*((vf - xl)*df - f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) {
			dxold = dx;
			dx = 0.5*(xh - xl);
			vf = xl + dx;
			if (xl == vf) return vf; //Change in root is negligible.
		}
		//Newton step acceptable. Take it.
		else {
			dxold = dx;
			dx = f / df;
			temp = vf;
			vf -= dx;
			if (temp == vf) return vf;
		}
		if (fabs(dx) < xacc) return vf; //Convergence criterion.

		//Recalculate the Rachford-Rice equation values. Pvap does not exist (=-1) for reaction families, which are species as defined currently
		f = 0.0;
		df = 0.0;
		//for (j = 0; j< NumMolecules; j++) {
		//	f += (Ith(y, j + 1) / sum) * ((Pvap[j] / pressure) - 1) / (1 + vf * ((Pvap[j] / pressure) - 1));
		//	//if (Pvap[j]>0.0) df += -(Ith(y, j + 1) / sum) * pow((Pvap[j] / gP) - 1, 2) / pow(1 + vf*((Pvap[j] / gP) - 1), 2);
		//	df += -(Ith(y, j + 1) / sum) * ((Pvap[j] / pressure) - 1)*((Pvap[j] / pressure) - 1) / ((1 + vf * ((Pvap[j] / pressure) - 1))*(1 + vf * ((Pvap[j] / pressure) - 1)));
		//}
		for (CI p = Components.begin(); p != Components.end(); ++p)
		{
			Component* comp = *p;
			int j = comp->m_nIndex;

			f += (Ith(y, j + 1) / sum) * ((Pvap[j] / pressure) - 1) / (1 + vf * ((Pvap[j] / pressure) - 1));
			//if (Pvap[j]>0.0) df += -(Ith(y, j + 1) / sum) * pow((Pvap[j] / gP) - 1, 2) / pow(1 + vf*((Pvap[j] / gP) - 1), 2);
			df += -(Ith(y, j + 1) / sum) * ((Pvap[j] / pressure) - 1)*((Pvap[j] / pressure) - 1) / ((1 + vf * ((Pvap[j] / pressure) - 1))*(1 + vf * ((Pvap[j] / pressure) - 1)));
		}

		//The one new function evaluation per iteration.
		if (f < 0.0) //Maintain the bracket on the root.
			xl = vf;
		else
			xh = vf;
	}

	return -1.0; //Never get here, hopefully. This means no solution found. Practically should converge in very few (2-5) iterations.
}

