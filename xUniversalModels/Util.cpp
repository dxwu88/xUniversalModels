#include "stdafx.h"
#include "Constants.h"
#include "Component.h"
#include "ComponentCollection.h"
#include "Util.h"

Util::Util()
{
}

Util::~Util()
{
}

//!This subroutine does an atomic balance on the prod array
//!based on the Feed array.
//
//use CompData, only: CR, HY, i_H2            !Carbon, hydrogen array, hydrogen component index
bool Util::AtomBalance(double feedStreamMoles[], double prodStreamMoles[], double& CError, double &HError)
{
	// Calculate the atomic species in the feed
	double FeedC = 0.0;
	double FeedH = 0.0;

	typedef std::vector<Component*> ::iterator CI;
	std::vector<Component*> Components = ComponentCollection::m_vecComponents;
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		FeedC += feedStreamMoles[i] * comp->m_fCR;
		FeedH += feedStreamMoles[i] * comp->m_fHY;
	}

	// Check if no carbon in feed
	if (FeedC <= 0.0) return false;

	// Force the atomic balance
	double ProdC = 0.0;
	double ProdH = 0.0;
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		ProdC += prodStreamMoles[i] * comp->m_fCR;
		ProdH += prodStreamMoles[i] * comp->m_fHY;
	}

	CError = FeedC - ProdC;

	// Check if no carbon in prod
	if (ProdC <= 0.0) return false;

	// Normalize carbon
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		prodStreamMoles[i] *= (FeedC / ProdC);
	}

	// Normalize hydrogen and Adjust H2
	HError = FeedH - ProdH;
	int H2Index = ComponentCollection::GetComponentIndex("Hydrogen");
	prodStreamMoles[H2Index] = prodStreamMoles[H2Index] + (HError / 2.0);

	// Check for less than 0.0 and fix
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		if (prodStreamMoles[i] < 0.0) prodStreamMoles[i] = 0.0;
	}
}

double Util::TestFunction(double x)
{
	return (x*x*x - 4 * x - 9);
}

// declare a function: auto f = [](double x){ return (x*x*x -4*x - 9); };
//
double Util::brents_fun(std::function<double(double)> f, double lower_bound, double upper_bound, double TOL, int MAX_ITER)
{
	double a = lower_bound;
	double b = upper_bound;
	double fa = f(a);   // calculated now to save function calls
	double fb = f(b);   // calculated now to save function calls
	double fs = 0;      // initialize 

	// Signs of f(lower_bound) and f(upper_bound) must be opposites, throws exception if root isn't bracketed
	if (!(fa * fb < 0.0))
	{
		return DBL_MAX;
	}

	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
	if (abs(fa) < abs(b)) 
	{
		std::swap(a, b);
		std::swap(fa, fb);
	}

	double c = a;           // c now equals the largest magnitude of the lower and upper bounds
	double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
	bool mflag = true;      // boolean flag used to evaluate if statement later on
	double s = 0;           // Our Root that will be returned
	double d = 0;           // Only used if mflag is unset (mflag == false)

	for (int iter = 0; iter < MAX_ITER; ++iter)
	{
		// stop if converged on root or error is less than tolerance
		if (abs(b - a) < TOL)
		{
			//std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
			return s;
		} // end if

		if (fa != fc && fb != fc)
		{
			// use inverse quadratic interopolation
			s = (a * fb * fc / ((fa - fb) * (fa - fc))) 
				+ (b * fa * fc / ((fb - fa) * (fb - fc)))
				+ (c * fa * fb / ((fc - fa) * (fc - fb)));
		}
		else
		{
			// secant method
			s = b - fb * (b - a) / (fb - fa);
		}

		/*
		Crazy condition statement!:
		-------------------------------------------------------
		(condition 1) s is not between  (3a+b)/4  and b or
		(condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
		(condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
		(condition 4) (mflag is set and |b−c| < |TOL|) or
		(condition 5) (mflag is false and |c−d| < |TOL|)
		*/
		if (((s < (3 * a + b) * 0.25) || (s > b)) || (mflag && (abs(s - b) >= (abs(b - c) * 0.5))) 
			|| (!mflag && (abs(s - b) >= (abs(c - d) * 0.5))) || (mflag && (abs(b - c) < TOL)) 
			|| (!mflag && (abs(c - d) < TOL)))
		{
			// bisection method
			s = 0.5 * (a + b);
			mflag = true;
		}
		else
		{
			mflag = false;
		}

		fs = f(s);  // calculate fs
		d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
		c = b;      // set c equal to upper bound
		fc = fb;    // set f(c) = f(b)

		if (fa * fs < 0)   // fa and fs have opposite signs
		{
			b = s;
			fb = fs;    // set f(b) = f(s)
		}
		else
		{
			a = s;
			fa = fs;    // set f(a) = f(s)
		}

		// if magnitude of fa is less than magnitude of fb
		if (abs(fa) < abs(fb)) 
		{
			std::swap(a, b);     // swap a and b
			std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
		}
	} 

	// The solution does not converge or iterations are not sufficient
	return DBL_MAX;
} 
