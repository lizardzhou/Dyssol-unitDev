#define DLL_EXPORT
#include "Unit.h"
#include "math.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// Spray nozzle with 2 different models: 1) simple pressure; 2) pneumatic
/// Calculates Sauter diameter of outlet droplets and plot their size distribution (log-normal distribution)

CUnit::CUnit()
{
	// Basic unit's info
	m_sUnitName = "Spray nozzle";
	m_sAuthorName = "XYZhou";
	m_sUniqueID = "80D25FB2-6603-4EF2-BB81-F33F25DCEAE4";

	// Add ports
	AddPort("InPort", INPUT_PORT);
	AddPort("OutPort", OUTPUT_PORT);

	// Add groups and unit parameters
	AddGroupParameter("Model", simplePressure, { simplePressure, pneumatic }, { "Simple pressure nozzle", "Pneumatic nozzle" }, "Atomization model");
	AddConstParameter("D", 0, 2, 1, "mm", "Nozzle diameter"); // nozzle diameter
	AddConstParameter("deltaP", 30, 200, 50, "bar", "Pressure difference of nozzle"); // pressure drop in nozzle
	AddConstParameter("geoStdDev", 1, 100, 2.5, "-", "Geometric standard deviation of droplet size distribution"); // geometric standard deviation
	AddParametersToGroup("Model", "Simple pressure nozzle", { "D", "deltaP", "geoStdDev" });
	AddParametersToGroup("Model", "Pneumatic nozzle", { "D", "deltaP", "geoStdDev" });
}

CUnit::~CUnit()
{

}

std::string CUnit::GetCompoundKey(std::string _sCompoundName)

{
	// Check if compound is defined in material database
	if (!IsCompoundNameDefined(_sCompoundName)) {
		RaiseError("Compound \"" + _sCompoundName + "\" is not available in material database.");
	}

	// Return the specific compound key. Input value is the name of this compound in string type
	std::vector<std::string> vCompoundNames = GetCompoundsNames();
	for (unsigned i = 0; i < GetCompoundsNumber(); i++) {
		if (vCompoundNames[i] == _sCompoundName) {
			return GetCompoundsList().at(i);
		}
	}
	return "";
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

	// Get selected model
	m_model = static_cast<EModels>(GetGroupParameterValue("Model"));

	//Initialization for different models
	switch (m_model) {
		case simplePressure:	InitializeSimplePressure(_dTime); 
								break;
		case pneumatic:			InitializePneumatic(_dTime); 
								break;
	}
}

void CUnit::InitializeSimplePressure(double _dTime) {
	//Add plot for droplet size distribution
	AddPlot("Plot 1", "Droplet size [mm]", "PDF");
	AddCurveOnPlot("Plot 1", "Probability density function");
	AddPlot("Plot 2", "Droplet size [mm]", "CDF");
	AddCurveOnPlot("Plot 2", "Cumulative distribution function");

	//Add results: calculated Sauter, median and mean diameter
	AddStateVariable("Sauter diameter", 0, true);
	AddStateVariable("Count median diameter", 0, true);
	AddStateVariable("Count mean diameter", 0, true);
	AddStateVariable("Count mode diameter", 0, true);
}

void CUnit::InitializePneumatic(double _dTime) {
	//Add plot for droplet size distribution
	AddPlot("Plot 1", "Droplet size [mm]", "PDF");
	AddCurveOnPlot("Plot 1", "Probability density function");
	AddPlot("Plot 2", "Droplet size [mm]", "CDF");
	AddCurveOnPlot("Plot 2", "Cumulative distribution function");

	//Add results: calculated Sauter, median and mean diameter
	AddStateVariable("Sauter diameter", 0, true);
	AddStateVariable("Count median diameter", 0, true);
	AddStateVariable("Count mean diameter", 0, true);
	AddStateVariable("Count mode diameter", 0, true);
}

void CUnit::Simulate(double _dTime)
{
	//Simulation for different models
	switch (m_model) {
		case simplePressure:	SimulateSimplePressure(_dTime);
								break;
		case pneumatic:			SimulatePneumatic(_dTime);
								break;
	}

	//Save added state variables for the output
	SaveStateVariables(_dTime);
}

void CUnit::SimulateSimplePressure(double _dTime) {

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
	double massLiquid = pInStream->GetMassFlow(_dTime);
	double densityLiquid = pInStream->GetCompoundTPDProp(_dTime, compounds[0], DENSITY);
	double volumeLiquid = massLiquid / densityLiquid;
	double sigmaLiquid = GetCompoundConstant(compounds[0], CONST_PROP_USER_DEFINED_01);
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, compounds[0], VISCOSITY);

	//Physical properties of inlet gas
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, compounds[1], DENSITY);

	//Unit parameters
	double D = GetConstParameterValue("D") * 1e-3; // Convert unit [mm] to [m]
	double deltaP = GetConstParameterValue("deltaP") * 1e5; // Convert unit [bar] to [Pa]
	double gDev = GetConstParameterValue("geoStdDev");

	//Calculate Sauter-diameter of droplets in [mm]
	double sauter = 1e3 * 2.3 * D * pow((4 * volumeLiquid) / (pow(D, 2) * MATH_PI * sqrt(2 * deltaP / densityLiquid)), 0.25) *
									pow(deltaP * D / sigmaLiquid, -0.25) *
									pow(D * sqrt(deltaP * densityLiquid) / viscosityLiquid, -0.25) *
									pow(densityLiquid / densityGas, 0.25);
	//Calculate count median diameter in [mm]
	double median = sauter * exp(-2.5 * pow(log(gDev),2));
	//Calculate count mean diameter in [mm]
	double mean = median * exp(0.5 * pow(log(gDev), 2));
	//Calculate count mode diameter in [mm]
	double mode = median * exp(-0.5 * pow(log(gDev), 2));

	//Set calculated values for the output
	SetStateVariable("Sauter diameter", sauter);
	SetStateVariable("Count median diameter", median);
	SetStateVariable("Count mean diameter", mean);
	SetStateVariable("Count mode diameter", mode);
	
	//Plot PDF
	std::vector<double> diameter1;
	std::vector<double> pdf;
	for (int i = 0; i < 200; i++) {
		diameter1.push_back(i*1e-3+1e-5);
		pdf.push_back((1 / (sqrt(2*MATH_PI)*log(gDev)*diameter1[i])) * exp(-0.5 * pow(log(diameter1[i]/median)/log(gDev),2)));
	}
	AddPointOnCurve("Plot 1", "Probability density function", diameter1, pdf);

	//Plot CDF
	std::vector<double> diameter2;
	std::vector<double> cdf;
	for (int i = 0; i < 200; i++) {
		diameter2.push_back(i*1e-3 + 1e-5);
		cdf.push_back(0.5 + 0.5 * erf((log(diameter2[i]/median))/(sqrt(2)*log(gDev))));
	}
	AddPointOnCurve("Plot 2", "Cumulative distribution function", diameter2, cdf);

	//Set outlet flow mass
	pOutStream->SetMassFlow(_dTime, massLiquid);
}

void CUnit::SimulatePneumatic(double _dTime) {
	CMaterialStream* pInStream = GetPortStream("InPort");
	CMaterialStream* pOutStream = GetPortStream("OutPort");

	//pOutStream->CopyFromStream(pInStream, _dTime);

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
	double sigmaLiquid = GetCompoundConstant(compounds[0], CONST_PROP_USER_DEFINED_01);
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, compounds[0], VISCOSITY);

	//Physical properties of inlet gas
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, compounds[1], DENSITY);

	//Unit parameters
	double D = GetConstParameterValue("D") * 1e-3; // Convert unit [mm] to [m]
	double deltaP = GetConstParameterValue("deltaP") * 1e5; // Convert unit [bar] to [Pa]
	double stdDev = GetConstParameterValue("stdDev") * 1e-3; // Convert unit [mm] to [m]
	double gDev = GetConstParameterValue("gDev");

	//Calculate Sauter-diameter of droplets in [mm]
	double sauter = 1e3 * 0.35 * D * pow(deltaP * D / (sigmaLiquid * pow(1 + massGas / massLiquid, 2)), -0.4) *
									 (1 + 2.5 * viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D));
	//Calculate count median diameter in [mm]
	double median = sauter * exp(-2.5 * pow(log(gDev), 2));
	//Calculate mean diameter in [mm]
	double mean = median * exp(0.5 * pow(stdDev, 2));

	//Plotting PDF
	

	//Plotting CDF
	
	
	SetStateVariable("Sauter diameter", sauter);
	SetStateVariable("Count median diameter", median);
	SetStateVariable("Mean diameter", mean);

	pOutStream->SetMassFlow(_dTime, massLiquid + massGas);
}

void CUnit::Finalize()
{

}