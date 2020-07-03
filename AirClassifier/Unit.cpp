#define DLL_EXPORT
#include "Unit.h"
#include "math.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// An air classifiying process separating particles according to their sinking velocity in a fluid stream.

CUnit::CUnit()
{
	/// Basic unit's info ///
	m_sUnitName = "Air classifier";
	m_sAuthorName = "XYZhou";
	m_sUniqueID = "211D0E54C80A4F3EB464671EEA222932"; // DO NOT change the ID, because the unit parameter for simulation is connected to it!

	/// Add ports ///
	AddPort("Input", INPUT_PORT);
	AddPort("Coarse", OUTPUT_PORT);
	AddPort("Fines", OUTPUT_PORT);

	/// Add unit parameters ///
	AddConstParameter("A", 0.01, 100, 0, "Cross-setional area [m2]", "");		// A: cross-sectional area in m2
	AddConstParameter("x", 0.01, 10, 0, "Sharpness factor [-]", "");			// x: sharpness factor


	/// Add user data to model ///
	m_NLModel.SetUserData(this);
}

CUnit::~CUnit()
{

}

void CUnit::Initialize(double _dTime)
{
	/// Check Simulation Setup
	// Check for gas and solid phases
	if (!IsPhaseDefined(SOA_VAPOR))
		RaiseError("Gas phase is not defined! Simulation aborted.");
	if (!IsPhaseDefined(SOA_SOLID))
		RaiseError("Solid phase is not defined! Simulation aborted.");

	// Check for size distribution
	if (!IsDistributionDefined(DISTR_SIZE))
		RaiseError("Particle size distribution not defined! Simulation aborted.");

	/// Clear all state variables in model ///
	m_NLModel.ClearVariables();

	/// Get number of diameter classes
	unsigned num_classes = GetClassesNumber(DISTR_SIZE);

	/// Add variable to the model of nonlinear equation system for each particle size ///
	for (unsigned i = 0; i < num_classes; ++i)
		m_NLModel.AddNLVariable(1.0);		// v_rel_i (relative velocity for each particle size class)

	/// Set model to the solver ///
	if (!m_NLSolver.SetModel(&m_NLModel))
		RaiseError(m_NLSolver.GetError());

	/// Add Plot
	AddPlot("Plot", "Diameter [m]", "Separation efficiency [-]", "Time [s]");

}

void CUnit::Simulate(double _dTime)
{
	/// Set model's current simulation time ///
	m_NLModel.time = _dTime;

	/// Run solver ///
	if (!m_NLSolver.Calculate(_dTime, KIN_FP)) // Here the fixed point iteration should be used and explicitly declared
		RaiseError(m_NLSolver.GetError());
}

void CUnit::Finalize()
{

}

void CUnit::SaveState()
{
	/// Save solver's state ///
	m_NLSolver.SaveState();
}

void CUnit::LoadState()
{
	/// Load solver's state ///
	m_NLSolver.LoadState();
}

//////////////////////////////////////////////////////////////////////////
/// Non-linear solver

void CMyNLModel::CalculateFunctions(double* _pVars, double* _pFunc, void* _pUserData)
{
	/// Get pointer to air classifier unit ///
	auto unit = static_cast<CUnit*>(_pUserData);

	/// Get pointers to streams ///
	CMaterialStream* inStream = unit->GetPortStream("Input");
	CMaterialStream* outStreamC = unit->GetPortStream("Coarse");
	CMaterialStream* outStreamF = unit->GetPortStream("Fines");

	/// Overall parameter ///
	double g = 9.81;			// graviational acceleration

	/// Get diameter classes and their number ///
	unsigned num_classes = unit->GetClassesNumber(DISTR_SIZE);
	std::vector<double> d = unit->GetClassesMeans(DISTR_SIZE);

	/// Get stream parameters ///
	double rho_solid = inStream->GetPhaseTPDProp(time, DENSITY, SOA_SOLID);
	double rho_gas = inStream->GetPhaseTPDProp(time, DENSITY, SOA_VAPOR);
	double eta_gas = inStream->GetPhaseTPDProp(time, VISCOSITY, SOA_VAPOR);

	/// Get value of variables (v_rel_i) at current iteration of solver ///
	std::vector<double> v_rel;
	for (unsigned i = 0; i < num_classes; ++i)
		v_rel.push_back(_pVars[i]); // insert value of _pVars from the the end of vector v_rel


	/// Calculation of new function values of relative velocity ///
	for (unsigned i = 0; i < num_classes; ++i)
	{
		/// Reynolds number of particle classes Re_i ///
		double Re_i = fabs(v_rel[i]) * d[i] * rho_gas / eta_gas;
		/// Drag coefficient of particle classes Cwp_i ///
		double Cwp_i = 24.0 / Re_i + 4.0 / sqrt(Re_i) + 0.4;
		/// Relative velocity ///
		double v_rel_update_i = sqrt((4.0 * rho_solid * d[i] * g) / (3.0 * rho_gas * Cwp_i));
		/// Update function value ///
		_pFunc[i] = v_rel_update_i;
	}

}

void CMyNLModel::ResultsHandler(double _dTime, double* _pVars, void* _pUserData)
{
	/// Get pointer to air classifier unit ///
	auto unit = static_cast<CUnit*>(_pUserData);

	/// Get pointers to streams ///
	CMaterialStream* inStream = unit->GetPortStream("Input");
	CMaterialStream* outStreamC = unit->GetPortStream("Coarse");
	CMaterialStream* outStreamF = unit->GetPortStream("Fines");

	/// Get diameter classes and their number ///
	std::vector<double> d = unit->GetClassesMeans(DISTR_SIZE);
	unsigned num_classes = unit->GetClassesNumber(DISTR_SIZE);

	/// Initialize output streams ///
	// Setting total mass flow to zero allows for only setting phase mass flows 
	// in the end of the unit (total mass flow will be calculated automatically)
	outStreamC->CopyFromStream(inStream, _dTime);
	outStreamC->SetMassFlow(_dTime, 0);
	outStreamF->CopyFromStream(inStream, _dTime);
	outStreamF->SetMassFlow(_dTime, 0);

	/// Setup transformation matrices ///
	CTransformMatrix TInputToCoarse(DISTR_SIZE, num_classes);
	CTransformMatrix TInputToFines(DISTR_SIZE, num_classes);

	/// Get parameters ///
	double A = unit->GetConstParameterValue("A");
	double x = unit->GetConstParameterValue("x");

	/// Get stream parameters ///
	double dm_solid = inStream->GetPhaseMassFlow(_dTime, SOA_SOLID);
	double rho_solid = inStream->GetPhaseTPDProp(_dTime, DENSITY, SOA_SOLID);
	double dm_gas = inStream->GetPhaseMassFlow(_dTime, SOA_VAPOR);
	double rho_gas = inStream->GetPhaseTPDProp(_dTime, DENSITY, SOA_VAPOR);
	std::vector<double> wIn = inStream->GetPSD(_dTime, PSD_MassFrac);

	/// Calculate cut velocity ///
	double w_cut = dm_gas / (rho_gas * A);

	/// Calculate separation efficiency ///
	// Fraction of mass in coarse stream
	double wC_acc = 0;
	// Separation efficiency for each particle class
	std::vector<double> xiC;
	for (unsigned i = 0; i < num_classes; ++i)
	{
		/// Get value of variables (v_rel_i) after convergence of solver ///
		double v_rel_i = _pVars[i];

		/// Temporary value for separation of particle class to coarse stream ///
		double  xiC_i;

		/// Check values of relative velocity ///
		// If v_rel_i < 0, particles are faster than fluid, i.e. they will go to Fines
		// Else calculate separation based on functions
		if (v_rel_i < 0)
		{
			xiC_i = 0; // go to Fines -> fraction in Coarse is 0
		}
		else
		{
			xiC_i = 1 / (1 + w_cut * exp(x * (1 - pow(v_rel_i / w_cut, 3))) / v_rel_i);
		}

		/// Update fraction of mass that goes to coarse stream ///
		wC_acc += wIn[i] * xiC_i;

		/// Update transformation matrices of the separation ///
		// add values of separation efficiency on TM diagonal
		TInputToCoarse.SetValue(i, i, xiC_i);
		TInputToFines.SetValue(i, i, 1 - xiC_i);

		/// Save temporary separation value to vector ///
		xiC.push_back(xiC_i);
	}

	/// Set properties coarse stream ///
	// Apply transformation matrix to coarse stream
	outStreamC->ApplyTM(_dTime, TInputToCoarse);
	// Set coarse solid mass flow
	outStreamC->SetPhaseMassFlow(_dTime, SOA_SOLID, dm_solid * wC_acc);

	/// Set properties fines stream ///
	// Apply tranformation matrix to fines stream
	outStreamF->ApplyTM(_dTime, TInputToFines);
	// Set gas mass flow
	outStreamF->SetPhaseMassFlow(_dTime, SOA_VAPOR, dm_gas);
	// Set solid mass flow
	outStreamF->SetPhaseMassFlow(_dTime, SOA_SOLID, dm_solid * (1 - wC_acc));

	/// Plotting
	// Separation efficiency (to coarse)
	unit->AddCurveOnPlot("Plot", _dTime);
	unit->AddPointOnCurve("Plot", _dTime, d, xiC);
}