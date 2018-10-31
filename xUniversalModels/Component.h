#pragma once
class Component
{
public:
	Component(int index, std::string strName);
	Component(int index, std::string strName, PhaseEnum phase, double MW, double SG,
		double idealGasEnthalpyFactors[], EnthalpyEquationTypeEnum eete, double idealGasGibbsFactors[],
		double ADS, double CR, double HY, double OX, double NI, double SU);

	~Component();

public:
	int m_nIndex;
	std::string m_strName;
	PhaseEnum m_ePhase;
	double m_fMW;
	double m_fSG;
	
	// C
	double m_fCR;
	// H
	double m_fHY;
	// O
	double m_fOX;
	// N
	double m_fNI;
	// S
	double m_fSU;

	EnthalpyEquationTypeEnum m_eEnthalpyEquationType;
	double m_fIdealgasEnthalpyFactors[7];

	// gibbs
	double m_fIdealGasGibbsFactors[5];

	double m_fADS;

};

