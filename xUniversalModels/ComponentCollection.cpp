#include "stdafx.h"
#include "ComponentCollection.h"

boost::unordered_map<int, Component*> ComponentCollection::m_mapComponents;
std::vector<Component*> ComponentCollection::m_vecComponents;

ComponentCollection::ComponentCollection()
{
	int index;
	std::string strName;
	PhaseEnum phase;
	double MW;
	double SG;
	double idealgasEnthalpyFactors[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	EnthalpyEquationTypeEnum eete;
	double idealGasGibbsFactors[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	double ADS;
	double CR;
	double HY;
	double OX;
	double NI;
	double SU;

	index = 0;
	strName = "Hydrogen";
	phase = GAS;
	MW = 2.015880108;
	SG = 69.85910034;
	double E[] = { 13.69972346, 4.742345521, 2466.0, 1.865190289, 567.6, 0.0, -29767.5 };
	memcpy_s(idealgasEnthalpyFactors, 7, E, 7);
	eete = DIPPR117;
	double G[] = { -32.67042075, 0.000381864, -6.90937E-07, 5.19185E-10, -1.38451E-13 };
	memcpy_s(idealGasGibbsFactors, 5, G, 5);
	ADS = 1.0;
	CR = 0;
	HY = 2;
	OX = 0;
	NI = 0;
	SU = 0;
	Component* comp = new Component(index, strName, phase, MW, SG, idealgasEnthalpyFactors, eete, 
									idealGasGibbsFactors, ADS, CR, HY, OX, NI, SU);
	m_vecComponents.push_back(comp);

	// add all components into the map
	typedef std::vector<Component*> ::iterator CI;
	for (CI p = m_vecComponents.begin(); p != m_vecComponents.end(); ++p)
	{
		Component* comp = *p;
		m_mapComponents[comp->m_nIndex] = comp;
	}
}

ComponentCollection::~ComponentCollection()
{
	// delete all components in vector
	typedef std::vector<Component*> ::iterator CI;
	for (CI p = m_vecComponents.begin(); p != m_vecComponents.end(); ++p)
	{
		Component* comp = *p;
		delete comp;
	}

	m_vecComponents.clear();
}

int ComponentCollection::GetComponentIndex(std::string compName)
{
	typedef std::vector<Component*> ::iterator CI;
	for (CI p = m_vecComponents.begin(); p != m_vecComponents.end(); ++p)
	{
		Component* comp = *p;
		if (comp->m_strName.compare(compName) == 0)
		{
			return comp->m_nIndex;
			break;
		}
	}
	return -1;
}
