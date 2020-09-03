#pragma once

#include "UnitDevelopmentDefines.h"

enum EModels : size_t {
	singleFluid, twoFluidExternal, twoFluidInternal, rotary
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
	void SimulateSingleFluid(double _dTime);
	void SimulateTwoFluidExternal(double _dTime);
	void SimulateTwoFluidInternal(double _dTime);
	void SimulateRotary(double _dTime);
	void Finalize() override;
};