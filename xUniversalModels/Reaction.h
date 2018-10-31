#pragma once
#include "Component.h"

class Reaction
{
public:
	Reaction();
	Reaction(int index, ReactorTypeEnum reactorType, ReactionTypeEnum reactionType, double rateConstant, double activeEnergy);
	~Reaction();

public:
	int m_nIndex;
	ReactorTypeEnum m_eReactorType;
	ReactionTypeEnum m_eReactionType;
	double m_fKEQRX;

	double m_fActiveEnergy;
	double m_fRateConstant;
	double m_fRate;

	// Component/Stoich pair vector
	std::vector<std::pair<Component*, double>> m_vecComponentStoich;
	std::vector<std::pair<Component*, double>> m_vecReactantComponentStoich;
	std::vector<std::pair<Component*, double>> m_vecProductComponentStoich;

public:
	double GetKEQRX() {	return m_fKEQRX; }
	void SetKEQRX(double KEQRX) { m_fKEQRX = KEQRX; }

};

