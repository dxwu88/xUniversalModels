#pragma once
#include "Reaction.h"

class ReactionCollection
{
public:
	ReactionCollection();
	~ReactionCollection();

public:
	static std::vector<Reaction*> m_vecForwardReactions;
	static std::vector<Reaction*> m_vecReverseReactions;

private:
	double m_fVaporFraction;

public:
	double GetVaporFraction() { return m_fVaporFraction; }
	void SetVaporFraction(double vf) { m_fVaporFraction = vf; }
	bool CheckStoichiometry();

};

