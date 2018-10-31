#include "stdafx.h"
#include "Component.h"

Component::Component(int index, std::string strName)
{
	m_nIndex = index;
	m_strName = strName;
}

Component::Component(int index, std::string strName, PhaseEnum phase, double MW, double SG,
	double idealGasEnthalpyFactors[], EnthalpyEquationTypeEnum eete, double idealGasGibbsFactors[],
	double ADS, double CR, double HY, double OX, double NI, double SU)
{
	m_nIndex = index;
	m_strName = strName;
	m_ePhase = phase;
	m_fMW = MW;
	m_fSG = SG;

	m_fCR = CR;
	m_fHY = HY;
	m_fOX = OX;
	m_fNI = NI;
	m_fSU = SU;

	m_eEnthalpyEquationType = eete;
	memcpy_s(m_fIdealgasEnthalpyFactors, 7, idealGasEnthalpyFactors, 7);

	// gibbs
	memcpy_s(m_fIdealGasGibbsFactors, 5, idealGasGibbsFactors, 5);

	m_fADS = ADS;
}

Component::~Component()
{
}
