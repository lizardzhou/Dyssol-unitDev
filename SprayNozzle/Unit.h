#pragma once

#include "UnitDevelopmentDefines.h"

enum EModels : size_t {
	simplePressure, pneumaticExternal, pneumaticInternal, rotary
};

class CUnit : public CSteadyStateUnit
{
public:
	CUnit();
	~CUnit();

	std::string GetCompoundKey(std::string _sCompoundName);
	void calcDiameter(double &sauter, double &gDev, double &median, double &mean, double &mode);
	EModels m_model;
	void Initialize(double _dTime) override;
	void InitializeOhnesorge(double _dTime);	
	void InitializeWeber(double _dTime);
	void Simulate(double _dTime) override;
	void SimulateSimplePressure(double _dTime);
	void SimulatePneumaticExternal(double _dTime);
	void SimulatePneumaticInternal(double _dTime);
	void SimulateRotary(double _dTime);
	void Finalize() override;
};