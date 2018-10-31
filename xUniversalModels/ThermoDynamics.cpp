#include "stdafx.h"
#include "Constants.h"
#include "Component.h"
#include "ComponentCollection.h"
#include "Reaction.h"
#include "ReactionCollection.h"
#include "ThermoDynamics.h"

ThermoDynamics::ThermoDynamics()
{
}

ThermoDynamics::~ThermoDynamics()
{
}

// This subroutine calculates the ideal gas enthalpy of the stream
// in units of KJ / HR.
// 
// Variables:
// streamMoles(TotalComps)	real * 8      kg - moles / hr of material
// temperature              real * 8      Temperature, K
// return Enthalpy          real * 8      Ideal gas enthalpy of stream, kJ / hr
//!
//!**********************************************************************
//use CompData, only: TotalComps, AAE, BE, CE, DE, EE, FE, HF25, MW, EOFF, HEquationType
//
//double ThermoDynamics::IGEnthalpy(FeedMoles, temperature, Enthalpy)
double ThermoDynamics::IdealGasEnthalpy(double streamMoles[], double temperature)
{
	typedef std::vector<Component*> ::iterator CI;
	std::vector<Component*> Components = ComponentCollection::m_vecComponents;

	// Calculate Ideal gas enthalpy(KJ / KG)
	double T2 = temperature * temperature;
	double T3 = temperature * T2;
	double T4 = temperature * T3;
	double T5 = temperature * T4;
	
	double Enthalpy = 0.0;
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;

		// Emass, KJ / KG
		double EMass = 0.0;
		if (comp->m_eEnthalpyEquationType == POLYS)
		{
			// poly5 form of ethalpy equation(non - DIPPR components)
			EMass = comp->m_fIdealgasEnthalpyFactors[0] + temperature * comp->m_fIdealgasEnthalpyFactors[1] + T2 * comp->m_fIdealgasEnthalpyFactors[2] + T3 * comp->m_fIdealgasEnthalpyFactors[3] + T4 * comp->m_fIdealgasEnthalpyFactors[4] + T5 * comp->m_fIdealgasEnthalpyFactors[5];
		}
		else if (comp->m_eEnthalpyEquationType == DIPPR117)
		{
			// DIPPR 117 form of enthalpy equation
			EMass = comp->m_fIdealgasEnthalpyFactors[0] * temperature + comp->m_fIdealgasEnthalpyFactors[1] * comp->m_fIdealgasEnthalpyFactors[2] * (cosh(comp->m_fIdealgasEnthalpyFactors[2] / temperature) / sinh(comp->m_fIdealgasEnthalpyFactors[2] / temperature))
				- comp->m_fIdealgasEnthalpyFactors[3] * comp->m_fIdealgasEnthalpyFactors[4] * tanh(comp->m_fIdealgasEnthalpyFactors[4] / temperature) + comp->m_fIdealgasEnthalpyFactors[5];
		}
		
		Enthalpy += streamMoles[comp->m_nIndex] * (comp->m_fMW * EMass + comp->m_fIdealgasEnthalpyFactors[6]);
	}

	return Enthalpy;
}

//!This function calculates Gibbs free energy of formation for component
//!"ICOMP".The units of the returned value are in kJ / kg - mole.
//!
//!Variables:
//!ICOMP    integer * 4   the component number
//!temperature        real * 8      Temperature, deg K
//!
//!Method :
//	!The calculation is done by using the curve fit polynomial equation
//	!in temperature stored in compdata.
//	!
//	!GIBBS = A + B * T + C * T*T + D * T*T*T + E * T*T*T*T
double ThermoDynamics::GIBBS(Component* comp, double temperature)
{
	// Calculate the GIBBS energy by the curve fit in t
	double GIBBS = comp->m_fIdealGasGibbsFactors[0] + temperature * (comp->m_fIdealGasGibbsFactors[1] + temperature * (comp->m_fIdealGasGibbsFactors[2]
		+ temperature * (comp->m_fIdealGasGibbsFactors[3] + comp->m_fIdealGasGibbsFactors[4] * temperature)));

	return GIBBS;
}

//!This function calculates the ideal gas heat capacity in
//!kJ / (kmol*K) at temperture temperature.The DIPPR enthalpy equation
//!is used to estimate the heat capacity.
//!
//!Inputs
//!streamMoles(KL)       real * 8      Molar flow rate, Kg - mol / hr
//!temperature       real * 8      Temperture, K
//!
//!Return value
//!HeatCapacity    real * 8  Constant Pressure Heat Capacity, KJ / (kmol*K)
//!
//use CompData, only: TotalComps, AAE, BE, CE, DE, EE, FE, MW, HEquationType
double ThermoDynamics::HeatCapacity(double streamMoles[], double temperature)
{
	typedef std::vector<Component*> ::iterator CI;
	std::vector<Component*> Components = ComponentCollection::m_vecComponents;
	
	double T2 = temperature * temperature;
	double T3 = T2 * temperature;
	double T4 = T3 * temperature;

	double HeatCapacity = 0.0;
	double YTotal = 0.0;;
	for (CI p = Components.begin(); p != Components.end(); ++p)
	{
		Component* comp = *p;
		int i = comp->m_nIndex;

		if (streamMoles[i] > 0.0)
		{
			double CPMass = 0.0;
			YTotal += streamMoles[i];

			// Calculate Heat capacity on mass basis in KJ using mass enthalpy equation derivative w.r.t.T
			if (comp->m_eEnthalpyEquationType == POLYS)
			{
				CPMass = comp->m_fIdealgasEnthalpyFactors[1] + 2.0 * comp->m_fIdealgasEnthalpyFactors[2] *temperature + 3.0*comp->m_fIdealgasEnthalpyFactors[3] *T2
					+ 4.0*comp->m_fIdealgasEnthalpyFactors[4] *T3 + 5.0*comp->m_fIdealgasEnthalpyFactors[5] *T4;
			}
			else if (comp->m_eEnthalpyEquationType == DIPPR117)
			{
				CPMass = comp->m_fIdealgasEnthalpyFactors[0] + comp->m_fIdealgasEnthalpyFactors[1] * pow(comp->m_fIdealgasEnthalpyFactors[2] / temperature / sinh(comp->m_fIdealgasEnthalpyFactors[2] / temperature), 2.0)
					+ comp->m_fIdealgasEnthalpyFactors[3] * pow(comp->m_fIdealgasEnthalpyFactors[4] / temperature / cosh(comp->m_fIdealgasEnthalpyFactors[4] / temperature), 2.0);
			}

			HeatCapacity += streamMoles[i] * comp->m_fMW * CPMass;
		}
	}

	return HeatCapacity / YTotal;
}

//!This function calculates the reaction equilibrium constant at temp = t
//!for reaction no ''norxn''.The equilibrium constant will be used with molar
//!concentrations rather than partial pressures or mole fractions so an
//!adjustment is calculated.This is for vapor phase only.
//!
//!VARIABLES:
//!
//!NORXN    INT * 4       the reaction number
//!T        real * 8      Reaction temperature, deg K
//!PA       real * 8      Reaction pressure, atm
//!
//!METHOD:
//!
//!The following equation relates Gibbs energies and chemical eq.
//!
//!KEQ = EXP(-(DELTA GIBBS) / (RGasJ*T))
//!
//!where: GIBBS(i, T) = Gibbs free energy of component i at
//!temperature = T.
//!UNITS = KJ / KGMOLE
//!Gibbs is calculated for the vapor phase
//!
//!RGasJ = Gas constant, 8.314 kJ / (kgmole*deg K)
//!
//!For real gases the equilibrium constant is a ratio of fugacity rather
//!than pressure.The compressibility factor Z is used to convert KEQ to
//!a mole fraction basis.KEQ = K - FUGACITY / K - Z
//!
//!The equilibrium constant, KEQ, is defined for partial pressure of ideal
//!gases or fugacity.If it will be used with molar concentrations then it
//!needs to be corrected to this basis.
//!
//!yConc(i) = Mole fraction(i) * PA / (RGasM3 * temperature)
//!
//!KEQ - CONC = KEQ * (RGasM3*temperature)**(-Delta stoichiometric constants)
void ThermoDynamics::KEQRX(double temperature, double pressure)
{
	typedef std::vector<Reaction*> ::iterator CR;
	std::vector<Reaction*> forwardReactions = ReactionCollection::m_vecForwardReactions;

	typedef std::vector<std::pair<Component*, double>> ::iterator CS;

	// Calculate the delta Gibbs energy from feeds and products
	//double RGasJ = 8.3144621;
	double DELTAG = 0.0;
	double CFACT = 0.0;

	for (CR pr = forwardReactions.begin(); pr != forwardReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;
		
		std::vector<std::pair<Component*, double>> ComponentStoichs = reaction->m_vecComponentStoich;

		for (CS pcs = ComponentStoichs.begin(); pcs != ComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;
			
			double gibbs = GIBBS(pair.first, temperature);
			DELTAG += pair.second * gibbs;
			CFACT -= pair.second;
		}

		// Calculate the equilibrium constant from Gibbs free energy
		double KEQRX = exp(-(DELTAG) / (temperature * RGasJ));

		// Change for use with concentrations rather than mole fractions or fugacity
		KEQRX *= pow(RGasM3 * temperature, CFACT);
		reaction->SetKEQRX(KEQRX);
	}

	std::vector<Reaction*> reverseReactions = ReactionCollection::m_vecReverseReactions;

	for (CR pr = reverseReactions.begin(); pr != reverseReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;

		std::vector<std::pair<Component*, double>> ComponentStoichs = reaction->m_vecComponentStoich;

		for (CS pcs = ComponentStoichs.begin(); pcs != ComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;

			double gibbs = GIBBS(pair.first, temperature);
			DELTAG += pair.second * gibbs;
			CFACT -= pair.second;
		}

		// Calculate the equilibrium constant from Gibbs free energy
		double KEQRX = exp(-(DELTAG) / (temperature * RGasJ));

		// Change for use with concentrations rather than mole fractions or fugacity
		KEQRX *= pow(RGasM3 * temperature, CFACT);
		reaction->SetKEQRX(KEQRX);
	}
}

// single reaction KEQRX
double ThermoDynamics::KEQRX(Reaction* reaction, double temperature, double pressure)
{
	typedef std::vector<std::pair<Component*, double>> ::iterator CS;

	// Calculate the delta Gibbs energy from feeds and products
	//double RGasJ = 8.3144621;
	double DELTAG = 0.0;
	double CFACT = 0.0;

	std::vector<std::pair<Component*, double>> ComponentStoichs = reaction->m_vecComponentStoich;

	for (CS pcs = ComponentStoichs.begin(); pcs != ComponentStoichs.end(); ++pcs)
	{
		std::pair<Component*, double> pair = *pcs;

		DELTAG += pair.second * GIBBS(pair.first, temperature);
		CFACT -= pair.second;
	}

	// Calculate the equilibrium constant from Gibbs free energy
	double KEQRX = exp(-(DELTAG) / (temperature * RGasJ));

	// Change for use with concentrations rather than mole fractions or fugacity
	KEQRX *= pow(RGasM3 * temperature, CFACT);
	reaction->SetKEQRX(KEQRX);
	return KEQRX;
}
