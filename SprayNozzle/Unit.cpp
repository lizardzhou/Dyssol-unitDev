#define DLL_EXPORT
#include "Unit.h"
#include "math.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// Pneumatic spray nozzle
/// Calculates Sauter diameter of outlet droplets and plot their size distribution (normal distribution)

CUnit::CUnit()
{
	// Basic unit's info
	m_sUnitName = "Spray nozzle - Pneumatic";
	m_sAuthorName = "XYZhou";
	m_sUniqueID = "0E3F1ABC-583B-45DB-B980-99958B3D7467";

	/// Add ports
	AddPort("InPort", INPUT_PORT);
	AddPort("OutPort", OUTPUT_PORT);

	/// Add unit parameters
	//AddGroupParameter("Model", simplePressure, { simplePressure, peumatic }, { "Simple pressure nozzle", "Pneumatic nozzle" }, "Atomization model");
	AddConstParameter("D", 0, 2, 1, "mm", "Nozzle diameter"); // nozzle diameter
	AddConstParameter("deltaP", 30, 200, 50, "bar", "Pressure difference of nozzle"); // pressure drop in nozzle
	AddConstParameter("stdDev", 0, 100, 1e-3, "mm", "Standard deviaiton of droplet size distribution"); // standard deviation

	/*
	/// Allocation of parameter to each elememt in the model
	AddParametersToGroup("Model", "Simple pressure nozzle", {"D", "deltaP", "stdDev"});
	AddParametersToGroup("Model", "Pneumatic nozzle", { "D", "deltaP", "stdDev" });
	*/
}

CUnit::~CUnit()
{

}

void CUnit::Initialize(double _dTime)
{
	//Check for liquid and gas phases
	if (!IsPhaseDefined(SOA_LIQUID)) {
		RaiseError("Please define the liquid phase which is necessary for the spray process.");
	}
	if (!IsPhaseDefined(SOA_VAPOR)) {
		RaiseError("Please define the gas phase which is necessary for the spray process.");
	}

	//Add plot for droplet size distribution
}

void CUnit::Simulate(double _dTime)
{
	CMaterialStream* pInStream = GetPortStream("InPort");
	CMaterialStream* pOutStream = GetPortStream("OutPort");

	pOutStream->CopyFromStream(pInStream, _dTime);

	//Check compounds: currently only 2-compound model of liquid/gas available
	std::vector<std::string> compounds = GetCompoundsList();
	double soaComp0 = GetCompoundConstant(compounds[0], SOA_AT_NORMAL_CONDITIONS);
	double soaComp1 = GetCompoundConstant(compounds[1], SOA_AT_NORMAL_CONDITIONS);

	if (compounds.size() > 2) {
		RaiseError("Only 2-compound model is available.");
	}

	if (soaComp0 == 0 || soaComp1 == 0) {
		RaiseError("Solid phase is not applicable.");
	}
	else if (soaComp0 == soaComp1) {
		RaiseError("Compounds of same phase is not applicable.");
	}
	else if (soaComp0 == 2 && soaComp1 == 1) { // Sequence of compounds for calculation: 1st is liquid and 2nd is gas. Rearrange the position of liquid and gas compound if needed
		std::string exchange = compounds[0];
		compounds[0] = compounds[1];
		compounds[1] = exchange;
	}
	
	//Physical properties of inlet liquid
	double massLiquid = pInStream->GetCompoundMassFlow(_dTime, compounds[0], SOA_LIQUID, BASIS_MASS);
	double massGas = pInStream->GetCompoundMassFlow(_dTime, compounds[1], SOA_LIQUID, BASIS_MASS);
	double densityLiquid = pInStream->GetCompoundTPDProp(_dTime, compounds[0], DENSITY);
	double volumeLiquid = massLiquid / densityLiquid;
	double sigma = pInStream->GetCompoundInteractionProp(_dTime, compounds[0], compounds[1], INTERFACE_TENSION); //Interphase tension liquid/gas
						 //other option: pInStream->GetCompoundTPDProp(_dTime, compounds[0], TP_PROP_USER_DEFINED_01); // Surface tension of liquid
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, compounds[0], VISCOSITY);

	//Physical properties of inlet gas
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, compounds[1], DENSITY);

	//Unit parameters
	double D = GetConstParameterValue("D") * 1e-3; // Convert unit of nozzle diameter to [m]
	double deltaP = GetConstParameterValue("deltaP") * 1.013e5; // Convert unit of pressure difference to [Pa]
	double stdDev = GetConstParameterValue("stdDev") * 1e-3; // Convert unit of standard deviation to [m]

	//Calculate Sauter-diameter of droplets
	double sauterFlowOut = 0.35 * pow(deltaP * D / (sigma * pow(1 + massGas / massLiquid, 2)), -0.4) * 
						   (1 + 2.5 * viscosityLiquid / sqrt(sigma * densityLiquid * D));
	
	//Plotting calculated Sauter-diameter with user-input standard deviation

}

void CUnit::Finalize()
{

}