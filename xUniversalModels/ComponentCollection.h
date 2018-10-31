#pragma once
#include "Component.h"

class ComponentCollection
{
public:
	ComponentCollection();
	~ComponentCollection();

public:
	static boost::unordered_map<int, Component*> m_mapComponents;
	static std::vector<Component*> m_vecComponents;

public:
	static int GetComponentIndex(std::string compName);
};

