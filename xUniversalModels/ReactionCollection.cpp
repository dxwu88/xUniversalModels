#include "stdafx.h"
#include "Component.h"
#include "ComponentCollection.h"
#include "Reaction.h"
#include "ReactionCollection.h"

std::vector<Reaction*> ReactionCollection::m_vecForwardReactions;
std::vector<Reaction*> ReactionCollection::m_vecReverseReactions;

ReactionCollection::ReactionCollection()
{
	int index;
	ReactorTypeEnum reactorType;
	ReactionTypeEnum reactionType;
	double activeEnergy;
	double rateConstant;

	// reaction #0
	index = 0;
	reactorType = PFR;
	reactionType = LHHW;
	activeEnergy = 21000.0;
	rateConstant = 1.5;
	Reaction* reaction = new Reaction(index, reactorType, reactionType, rateConstant, activeEnergy);

	std::pair<Component*, double> pair = std::pair<Component*, double>(ComponentCollection::m_vecComponents[0], -2.0);
	reaction->m_vecComponentStoich.push_back(pair);
	//pair = std::pair<Component*, double>(ComponentCollection::m_vecComponents[1], -1.0);
	//reaction->m_vecComponentStoich.push_back(pair);
	//pair = std::pair<Component*, double>(ComponentCollection::m_vecComponents[2], 1.0);
	//reaction->m_vecComponentStoich.push_back(pair);

	m_vecForwardReactions.push_back(reaction);

	// reaction #1
	// ......

}

ReactionCollection::~ReactionCollection()
{
}

bool ReactionCollection::CheckStoichiometry()
{
	typedef std::vector<Reaction*> ::iterator CR;
	std::vector<Reaction*> forwardReactions = ReactionCollection::m_vecForwardReactions;

	typedef std::vector<std::pair<Component*, double>> ::iterator CS;
	
	for (CR pr = forwardReactions.begin(); pr != forwardReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;

		double C1 = 0.0;
		double H1 = 0.0;
		double O1 = 0.0;
		double N1 = 0.0;
		double S1 = 0.0;

		std::vector<std::pair<Component*, double>> ComponentStoichs = reaction->m_vecComponentStoich;

		for (CS pcs = ComponentStoichs.begin(); pcs != ComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;
			Component* comp = pair.first;

			C1 += pair.second * comp->m_fCR;
			H1 += pair.second * comp->m_fHY;
			O1 += pair.second * comp->m_fOX;
			N1 += pair.second * comp->m_fNI;
			S1 += pair.second * comp->m_fSU;
		}

		if (C1 != 0.0 || H1 != 0.0 || O1 != 0.0 || N1 != 0.0 || S1 != 0.0)
			return false;
	}

	std::vector<Reaction*> reverseReactions = ReactionCollection::m_vecReverseReactions;

	for (CR pr = reverseReactions.begin(); pr != reverseReactions.end(); ++pr)
	{
		Reaction* reaction = *pr;

		double C1 = 0.0;
		double H1 = 0.0;
		double O1 = 0.0;
		double N1 = 0.0;
		double S1 = 0.0;

		std::vector<std::pair<Component*, double>> ComponentStoichs = reaction->m_vecComponentStoich;

		for (CS pcs = ComponentStoichs.begin(); pcs != ComponentStoichs.end(); ++pcs)
		{
			std::pair<Component*, double> pair = *pcs;
			Component* comp = pair.first;

			C1 += pair.second * comp->m_fCR;
			H1 += pair.second * comp->m_fHY;
			O1 += pair.second * comp->m_fOX;
			N1 += pair.second * comp->m_fNI;
			S1 += pair.second * comp->m_fSU;
		}

		if (C1 != 0.0 || H1 != 0.0 || O1 != 0.0 || N1 != 0.0 || S1 != 0.0)
			return false;
	}

	return true;
}

