#define DLL_EXPORT
#include "Unit.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// A unit testing the basic functionality of classes BaseUnit, MaterialStream and Holdup

CUnit::CUnit()
{
	/// Basic unit's info ///
	m_sUnitName = "Basics";
	m_sAuthorName = "XYZhou";
	m_sUniqueID = "906354E1-ABCB-4575-8CFA-D59B3090D5FA";

	/// Add ports ///
	AddPort("InPort", INPUT_PORT);
	AddPort("OutPort", OUTPUT_PORT);

	/// Add unit parameters ///
	AddTDParameter("ParamTD", 0, 1e+6, 0, "kg", "Unit parameter description");
	AddConstParameter("ParamConst", 0, 1e+6, 0, "s", "Unit parameter description");
	AddStringParameter("ParamString", "Initial value", "Unit parameter description");

	/// Add holdups ///
	AddHoldup("Holdup");
}

CUnit::~CUnit()
{

}

void CUnit::Initialize(double _dTime)
{
	/// Add state variables ///
	AddStateVariable("VarName", 0, true);

	// Warning of liquid or vapor phase is not defined
	if ((!IsPhaseDefined(SOA_LIQUID)) || (!IsPhaseDefined(SOA_VAPOR))) {
		RaiseWarning("The state of aggregation is not defined! Please define liquid or vapor phase.");
	}

	// Add an internal stream
	AddMaterialStream("bufStream");

	// Add plot
	AddPlot("Time dependence of holdup mass", "Time [s]", "Mass [kg]");
	AddCurveOnPlot("Time dependence of holdup mass", "Curve1");
}

void CUnit::Simulate(double _dStartTime, double _dEndTime)
{
	/// Get pointers to streams ///
	CMaterialStream* pInStream = GetPortStream("InPort");
	CMaterialStream* pOutStream = GetPortStream("OutPort");
	CMaterialStream* bufStream = GetMaterialStream("bufStream");

	/// Get pointers to holdups ///
	CHoldup* pHoldup = GetHoldup("Holdup");

	// Add start time point to bufStream
	bufStream->AddTimePoint(_dStartTime);

	// Copy inlet stream into bufStream
	bufStream->CopyFromStream(pInStream, _dEndTime);

	// Set mass flow 12.5 kg/s of liquid phase in bufStream at time point 10 s
	bufStream->SetPhaseMassFlow(10, SOA_LIQUID, 12.5, BASIS_MASS);

	// Add inlet to the holdup on entire time interval
	pHoldup->AddStream(pInStream, _dStartTime, _dEndTime);

	// Copy the holdup into outlet stream at end time point with mass flow 1 kg/s
	pOutStream->CopyFromHoldup(pHoldup, _dStartTime, 1); // copy the holdup info at the beginning which is deleted afterwards
	pOutStream->CopyFromHoldup(pHoldup, _dEndTime, 1); // copy the holdup infor at the end which is deleted afterwards

	//pOutStream->AddTimePoint(_dEndTime, _dStartTime);

	// Set new temperature 320 K to outlet at time point 15 s
	//pOutStream->AddTimePoint(15);
	pOutStream->SetTemperature(_dEndTime, 320);

	// Plot holdup mass for all defined time points
	std::vector<double> times = GetAllDefinedTimePoints(_dStartTime, _dEndTime);
	for (int i = 0; i < times.size(); i++) {
		double x = times[i];
		double y = pHoldup->GetMass(times[i], BASIS_MASS);
		AddPointOnCurve("Time dependence of holdup mass", "Curve1", x, y);
	}

	/// Data acquisition
	// Get unit parameters
	double TDParameter = GetTDParameterValue("ParamTD", 5);
	double ConstParameter = GetConstParameterValue("ParamConst");
	std::string StringParameter = GetStringParameterValue("ParamString");
	// Get common compound information
	std::vector<std::string> compounds = GetCompoundsList(); //only one compound in task6, so only one element in compounds array
	double molarMass = GetCompoundConstant(compounds[0], MOLAR_MASS);
	double critTemp = GetCompoundConstant(compounds[0], CRITICAL_TEMPERATURE);
	double density = GetCompoundTPDProp(compounds[0], DENSITY, 273, 1e5);
	// Get tolerance
	double absTol = GetAbsTolerance();
	double relTol = GetRelTolerance();
	// Get overall properties of streams and holdups
	double massFlow = pInStream->GetMassFlow(2, BASIS_MASS);
	double massHoldup = pHoldup->GetMass(5, BASIS_MASS);
	//double outTemp1 = pOutStream->GetTemperature(15);
	double outTemp2 = pOutStream->GetTemperature(_dEndTime);
	double molarMassHoldup = pHoldup->GetOverallProperty(1, MOLAR_MASS);
	// Get solid distribution information
	std::vector<double> PSD_b3 = pHoldup->GetPSD(50, PSD_Q3);
	std::vector<double> PSD_s3 = pHoldup->GetPSD(50, PSD_q3);
}

void CUnit::SaveState()
{

}

void CUnit::LoadState()
{

}

void CUnit::Finalize()
{

}