#define DLL_EXPORT
#include "Unit.h"
#include "math.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// Simple pressure spray nozzle
/// Calculates Sauter diameter of outlet droplets and plot their size distribution (normal distribution)

CUnit::CUnit()
{
	// Basic unit's info
	m_sUnitName = "Spray nozzle - Simple";
	m_sAuthorName = "XYZhou";
	m_sUniqueID = "80D25FB2-6603-4EF2-BB81-F33F25DCEAE4";

	/// Add ports
	AddPort("InPort", INPUT_PORT);
	AddPort("OutPort", OUTPUT_PORT);

	/// Add unit parameters
	AddGroupParameter("Model", simplePressure, { simplePressure, peumatic }, { "Simple pressure nozzle", "Pneumatic nozzle" }, "Atomization model");
	AddConstParameter("D", 0, 2, 1, "mm", "Nozzle diameter"); // nozzle diameter
	AddConstParameter("deltaP", 30, 200, 50, "bar", "Pressure difference of nozzle"); // pressure drop in nozzle
	AddConstParameter("stdDev", 0, 100, 1e-3, "mm", "Standard deviaiton of droplet size distribution"); // standard deviation

	/// Allocation of parameter to each elememt in the model
	AddParametersToGroup("Model", "Simple pressure nozzle", {"D", "deltaP", "stdDev"});
	AddParametersToGroup("Model", "Pneumatic nozzle", { "D", "deltaP", "stdDev" });
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
	AddPlot("Plot 1", "Droplet size [mm]", "PDF");
	AddCurveOnPlot("Plot 1", "Curve 1");
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
	double densityLiquid = pInStream->GetCompoundTPDProp(_dTime, compounds[0], DENSITY);
	double volumeLiquid = massLiquid / densityLiquid;
	double sigmaLiquid = pInStream->GetCompoundInteractionProp(_dTime, compounds[0], compounds[1], INTERFACE_TENSION); //Interphase tension liquid/gas
						 //other option: pInStream->GetCompoundTPDProp(_dTime, compounds[0], TP_PROP_USER_DEFINED_01); // Surface tension of liquid
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, compounds[0], VISCOSITY);

	//Physical properties of inlet gas
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, compounds[1], DENSITY);

	//Unit parameters
	double D = GetConstParameterValue("D") * 1e-3; // Convert unit of nozzle diameter to [m]
	double deltaP = GetConstParameterValue("deltaP") * 1.013e5; // Convert unit of pressure difference to [Pa]
	double stdDev = GetConstParameterValue("stdDev");

	//Calculate Sauter-diameter of droplets
	double sauterFlowOut = 2.3 * ((4 * volumeLiquid) / (D * MATH_PI * sqrt(2 * deltaP / densityLiquid))) *
								 pow(D * sqrt(deltaP * densityLiquid / viscosityLiquid), -0.25) * 
								 pow(deltaP * D / sigmaLiquid, -0.25) *
								 pow(densityLiquid / densityGas, 0.25);
	
	//Plotting calculated Sauter-diameter with user-input standard deviation - now only a test
	std::vector<double> test = {1,2,3,4,5};
	std::vector<double> result;
	for (int i = 0; i < test.size(); i++) {
		result.push_back(pow(test[i],2));
	}
	AddPointOnCurve("Plot 1", "Curve 1", test, result);

}

void CUnit::Finalize()
{

}