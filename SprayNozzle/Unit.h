#pragma once

#include "UnitDevelopmentDefines.h"

/*
enum EModels : size_t {
	simplePressure, peumatic
};
*/

class CUnit : public CSteadyStateUnit
{
public:
	CUnit();
	~CUnit();

	void Initialize(double _dTime) override;
	void Simulate(double _dTime) override;
	void Finalize() override;
};