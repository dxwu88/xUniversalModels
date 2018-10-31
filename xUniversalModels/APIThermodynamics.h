#pragma once
#include "ThermoDynamics.h"
#include "Flash.h"

class APIThermodynamics : public ThermoDynamics, public Flash
{
public:
	APIThermodynamics();
	~APIThermodynamics();

	static double FunctionEnthalpy2Target(double temperature);
	int PHFlash(double streamMoles[], double InEnthalpy, double temperature, double PRESSatm, double StartT, double VF, int& Error);
	int CalcEnthalpy(double streamMoles[], double temperature, double PRatm, double Enthalpy, double VF);

};

