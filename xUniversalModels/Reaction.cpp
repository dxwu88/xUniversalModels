#include "stdafx.h"
#include "Reaction.h"

Reaction::Reaction()
{
	m_eReactorType = PFR;
	m_eReactionType = LHHW;
}

Reaction::Reaction(int index, ReactorTypeEnum reactorType, ReactionTypeEnum reactionType, double rateConstant, double activeEnergy)
{
	m_nIndex = index;

	m_eReactorType = reactorType;
	m_eReactionType = reactionType;
	m_fRateConstant = rateConstant;
	m_fActiveEnergy = activeEnergy;
}

Reaction::~Reaction()
{
}
