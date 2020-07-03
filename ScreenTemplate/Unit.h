#pragma once

#include "UnitDevelopmentDefines.h"

class CMyDAEModel : public CDAEModel
{
public:
	size_t m_nVariable0;

public:
	void CalculateResiduals(double _dTime, double* _pVars, double* _pDers, double* _pRes, void* _pUserData) override;
	void ResultsHandler(double _dTime, double* _pVars, double* _pDerivs, void* _pUserData) override;
};

class CUnit : public CDynamicUnit
{
private:
	CMyDAEModel m_Model;		// Model of DAE
	CDAESolver m_Solver;		// Solver of DAE

public:
	CUnit();
	~CUnit();

	void Initialize(double _dTime) override;
	void Simulate(double _dStartTime, double _dEndTime) override;
	void SaveState() override;
	void LoadState() override;
	void Finalize() override;
};