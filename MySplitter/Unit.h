#pragma once

#include "UnitDevelopmentDefines.h"

class CUnit : public CSteadyStateUnit
{
public:
	CUnit();
	~CUnit();

	void Initialize(double _dTime) override;
	void Simulate(double _dTime) override;
	void Finalize() override;
};