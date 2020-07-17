#pragma once

#include "UnitDevelopmentDefines.h"

enum EModels : size_t {
	simplePressure, pneumatic
};

class CUnit : public CSteadyStateUnit
{
public:
	CUnit();
	~CUnit();

	std::string GetCompoundKey(std::string _sCompoundName);
	EModels m_model;
	void Initialize(double _dTime) override;
	void InitializeSimplePressure(double _dTime);	
	void InitializePneumatic(double _dTime);
	void Simulate(double _dTime) override;
	void SimulateSimplePressure(double _dTime);
	void SimulatePneumatic(double _dTime);
	void Finalize() override;
};