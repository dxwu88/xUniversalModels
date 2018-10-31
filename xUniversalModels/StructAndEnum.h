#pragma once

enum ErrorCodeEnum
{
	SUCCESS = 0,
	Converge = 1,
	DIVIDEDBY0 = 2
};

enum PhaseEnum
{
	GAS = 0,
	LIQUID = 1,
	SOLID = 2,
};

enum FlashAlgorithmEnum
{
	EOS = 0,
	PR = 1
};

enum EnthalpyEquationTypeEnum
{
	POLYS = 0,
	DIPPR117 = 1
};

enum ReactorTypeEnum
{
	BATCH = 0,
	CSTR = 1,
	PFR = 2
};

enum ReactionTypeEnum
{
	LHHW = 0
};

enum TemperatureProfileEnum
{
	ISOTHERMAL = 0,
	ADIABATIC = 1,
	PROFILE = 2
};

typedef struct _StreamProperty
{
	int id;
} StreamPropertyStruct;

typedef struct _ProcessCondition
{
	int id;
} ProcessConditionStruct;

typedef struct _UserData
{
	int n;
	double x;

	double PressureAtm;
	double TempK;
	double TRef;
} UserDataStruct;

typedef struct _ReactorData
{
	int n;
	double x;

	double PressureAtm;
	double TempK;
	double TRef;

	TemperatureProfileEnum ReactorTemperatureProfile;
} ReactorDataStruct;