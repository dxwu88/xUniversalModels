#pragma once
#include "Component.h"
#include "Reaction.h"

class ThermoDynamics
{
public:
	ThermoDynamics();
	~ThermoDynamics();

public:
	static double IdealGasEnthalpy(double streamMoles[], double temperature);
	static double GIBBS(Component* comp, double temperature);
	static double HeatCapacity(double streamMoles[], double temperature);
	void KEQRX(double temperature, double pressure);
	static double KEQRX(Reaction* reaction, double temperature, double pressure);

};

