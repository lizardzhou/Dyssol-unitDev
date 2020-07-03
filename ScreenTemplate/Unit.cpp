#define DLL_EXPORT
#include "Unit.h"
#include "math.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// A dynamic screen model with a holdup. The separation sharpness reduces with time and also depends on the holdup mass.

CUnit::CUnit()
{
	/// Basic unit's info ///
	m_sUnitName = "Dynamic Screen";
	m_sAuthorName = "XYZhou";
	m_sUniqueID = "C7755DAF619C448D863D1CBCC13648BC"; // DO NOT change the ID, because the parameter of screen unit is connected to it!!!

	/// Add ports ///
	AddPort("Input", INPUT_PORT);
	AddPort("Coarse", OUTPUT_PORT);
	AddPort("Fines", OUTPUT_PORT);

	/// Add unit parameters ///
	AddConstParameter("alpha", 0, 100, 0, "Separation sharpness", "");	// alpha
	AddConstParameter("Xcut", 0, 1, 0, "Cut size", "");	 // Xcut
	AddConstParameter("Mout", 0, 100, 0, "Outlet mass flow [kg/s]", "");	// Mout
	AddConstParameter("k1", 0, 1, 0, "Time-dependent sharpness reduction factor [1/s]", "");	// k1
	AddConstParameter("k2", 0, 1, 0, "Time-dependent sharpness reduction factor [1/kg]", "");	// k2

	/// Add holdups ///
	AddHoldup("Holdup");

	/// Set this unit as user data of model ///
	m_Model.SetUserData(this);
}

CUnit::~CUnit()
{

}

void CUnit::Initialize(double _dTime)
{
	// Check flowsheet parameters
	if (!IsDistributionDefined(DISTR_SIZE)) {
		RaiseError("The size distribution is not defined! Simulation aborted.");
	}
	if (!IsPhaseDefined(SOA_SOLID)) {
		RaiseError("The solid phase is not defined! Simulation aborted.");
	}

	// Add plots
	AddPlot("Time dependence of separation sharpness", "Time [s]", "Sharpness [-]");
	AddCurveOnPlot("Time dependence of separation sharpness", "Curve");

	// Clear all state variables in model
	m_Model.ClearVariables();

	// Add state variables to a model 
	m_Model.AddDAEVariable(true, GetHoldup("Holdup")->GetMass(_dTime), 0);	// holdup mass: time-dependent
	m_Model.AddDAEVariable(true, GetConstParameterValue("alpha"), 0);	// separation sharpness: time-dependent
	m_Model.AddDAEVariable(false, GetConstParameterValue("Mout"), 0);	// outlet mass flow: steady

																		// Set tolerances to model ///
	m_Model.SetTolerance(GetRelTolerance() * 10, GetAbsTolerance() * 10);

	// Set model to a solver ///
	if (!m_Solver.SetModel(&m_Model))
		RaiseError(m_Solver.GetError());
}

void CUnit::Simulate(double _dStartTime, double _dEndTime)
{
	/// Run solver ///
	if (!m_Solver.Calculate(_dStartTime, _dEndTime))
		RaiseError(m_Solver.GetError());
}

void CUnit::SaveState()
{
	/// Save solver's state ///
	m_Solver.SaveState();
}

void CUnit::LoadState()
{
	/// Load solver's state ///
	m_Solver.LoadState();
}

void CUnit::Finalize()
{

}

//////////////////////////////////////////////////////////////////////////
/// DAE solver

void CMyDAEModel::CalculateResiduals(double _dTime, double* _pVars, double* _pDers, double* _pRes, void* _pUserData)
{
	// Get pointers to streams
	CUnit *unit = static_cast<CUnit*>(_pUserData);
	CMaterialStream *inStream = unit->GetPortStream("Input");	// Input
	CHoldup *holdup = unit->GetHoldup("Holdup");				// Holdup

																// Get time parameters
	double prevTime = holdup->GetLastTimePoint();
	double dTime = _dTime - prevTime; // Time difference

									  // Get values of input and internal parameters
	double k1 = unit->GetConstParameterValue("k1");		// k1
	double k2 = unit->GetConstParameterValue("k2");		// k2
	double mOut = unit->GetConstParameterValue("Mout");	// Mout
	double mIn = inStream->GetMassFlow(_dTime);			// Inlet mass flow at current time point
	double MhPrev = holdup->GetMass(prevTime);			// Holdup mass at previous time point

														// Calculate and set residuals
	double derMassHoldup = mIn - mOut;
	double derAlpha = -_pVars[1] * k1 - _pVars[1] * (MhPrev + derMassHoldup) * k2; // (MhPrev + derMassHoldup): holdup mass at current time point
	double valMassFlowOut; // valMassFlowOut = mCoarse + mFines;
	if (_pVars[0] > mOut*dTime) { // enough holdup mass: current holdup mass is larger than the amount flowing out in dTime period
		valMassFlowOut = mOut;
	}
	else // not enough holdup mass: current holdup mass is smaller than the amount flowing out in dTime period
	{
		valMassFlowOut = mIn;
	}
	_pRes[0] = _pDers[0] - derMassHoldup; // residual value of d(holdup mass)
	_pRes[1] = _pDers[1] - derAlpha; // residual value of d(alpha)
	_pRes[2] = _pVars[2] - valMassFlowOut; // residual value of outlet mass flow
}

void CMyDAEModel::ResultsHandler(double _dTime, double* _pVars, double* _pDerivs, void *_pUserData)
{
	// Get pointers to streams
	CUnit* unit = static_cast<CUnit*>(_pUserData);
	CMaterialStream* inStream = unit->GetPortStream("Input");		// Input
	CMaterialStream* outStreamC = unit->GetPortStream("Coarse");	// Coarse
	CMaterialStream* outStreamF = unit->GetPortStream("Fines");		// Fines
	CHoldup* holdup = unit->GetHoldup("Holdup");					// Holdup

	// Get values of unit parameters at current time point
	double xCut = unit->GetConstParameterValue("Xcut");
	double Mh = _pVars[0]; // holdup mass
	double alpha = _pVars[1]; // spearation sharpness
	double mFlowOut = _pVars[2]; // outlet mass flow

								 // Add points on plot
	unit->AddPointOnCurve("Time dependence of separation sharpness", "Curve", _dTime, alpha);

	// Mix input stream with holdup in the time interval [last time, current time]
	holdup->AddStream(inStream, holdup->GetLastTimePoint(), _dTime);

	// Obtain parameters for PSD calculation
	unsigned classesNum = unit->GetClassesNumber(DISTR_SIZE);	// number of grid classes for PSD
	std::vector<double> x = unit->GetClassesMeans(DISTR_SIZE);	// mean sizes of grid classes are listed in a vector with length classesNum
	std::vector<double> holdupPSD = holdup->GetPSD(_dTime, PSD_MassFrac);	// get holdup PSD as mass fractions at current time point: listed in a vector with length classesNum
	
	// Setup transformation matrices for one-dimensional size distribution
	CTransformMatrix THoldupToCoarse(DISTR_SIZE, classesNum);
	CTransformMatrix THoldupToFines(DISTR_SIZE, classesNum);

	// Calculate transformation matrices
	double massFracC = 0;
	for (unsigned i = 0; i < classesNum; i++) {
		for (unsigned j = 0; j < classesNum; j++) {
			if (i == j) {
				double value = 1 / (1 + pow(xCut / x[i], 2) * exp(alpha*(1 - pow(x[i] / xCut, 2))));
				THoldupToCoarse.SetValue(i, j, value);
				THoldupToFines.SetValue(i, j, 1 - value);
				massFracC += holdupPSD[i] * value;
			}
		}
	}

	// Copy holdup to output streams
	outStreamC->CopyFromHoldup(holdup, _dTime, mFlowOut * massFracC);
	outStreamF->CopyFromHoldup(holdup, _dTime, mFlowOut * (1 - massFracC));

	// Apply transformation matrix
	outStreamC->ApplyTM(_dTime, THoldupToCoarse);
	outStreamF->ApplyTM(_dTime, THoldupToFines);

	// Set new mass to the holdup
	holdup->SetMass(_dTime, Mh);
}