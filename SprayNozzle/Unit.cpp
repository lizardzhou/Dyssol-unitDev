#define DLL_EXPORT
#include "Unit.h"
#include "math.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// Spray nozzle with 3 different models: 1) simple pressure; 2) pneumatic with external mixing 3) pneumatic with internal mixing
/// User only needs to define the liquid phase. In case of pneumatic model, the gas phase will be calculated using ALR
/// Calculates Sauter diameter and count median, mean, mode diameter
/// Plot droplet size distribution (log-normal distribution)

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
	AddGroupParameter("Model", simplePressure, { simplePressure, pneumaticExternal, pneumaticInternal }, { "Simple pressure nozzle", "2-Phase nozzle with external mixing", "2-Phase nozzle with internal mixing" }, "Atomization model");
	AddConstParameter("D", 1, 2, 1, "mm", "Liquid oriface diameter"); // nozzle diameter for simple pressure model
	AddConstParameter("liquidD", 0.5, 2, 1, "mm", "Liquid oriface diameter"); // nozzle diameter for pneumatic pressure model
	AddConstParameter("ALR", 1e-5, 5, 1, "-", "Air-to-liquid mass flow ratio"); // air-to-liquid mass flow ratio for 2-phase nozzle
	AddConstParameter("liquidDeltaP", 30, 200, 50, "bar", "Pressure drop in nozzle"); //liquid pressure drop in simple pressure nozzle
	AddConstParameter("externalGasDeltaP", 0.5, 3, 3, "bar", "Pressure difference of nozzle"); // gas pressure drop in pneumatic nozzle with external mixing
	AddConstParameter("internalGasDeltaP", 0.5, 15, 3, "bar", "Pressure difference of nozzle"); // gas pressure drop in pneumatic nozzle with internal mixing
	AddConstParameter("geoStdDev", 1, 100, 2.5, "-", "Geometric standard deviation of droplet size distribution"); // geometric standard deviation
	AddStringParameter("Drop compound", "H2O", "Show the drop compound name");
	AddParametersToGroup("Model", "Simple pressure nozzle", { "D", "liquidDeltaP", "geoStdDev" });
	AddParametersToGroup("Model", "2-Phase nozzle with external mixing", { "liquidD", "externalGasDeltaP", "ALR", "geoStdDev" });
	AddParametersToGroup("Model", "2-Phase nozzle with internal mixing", { "liquidD", "internalGasDeltaP", "ALR", "geoStdDev" });
}

CUnit::~CUnit()
{

}

std::string CUnit::GetCompoundKey(std::string _sCompoundName){

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
	//Check for gas phase
	if (!IsPhaseDefined(SOA_VAPOR)) {
		RaiseError("Please define the gas phase which is necessary for the spray process.");
	}
	//check for solid phase applied for PSD calculation
	if (!IsPhaseDefined(SOA_SOLID)) {
		RaiseError("Please define the solid phase which is necessary for the particle size distribution.");
	}

	// Get selected model
	m_model = static_cast<EModels>(GetGroupParameterValue("Model"));

	//Initialization for different models
	switch (m_model) {
		case simplePressure:	InitializeSimplePressure(_dTime); 
								break;
		case pneumaticExternal:	InitializePneumaticExternal(_dTime); 
								break;
		case pneumaticInternal:	InitializePneumaticInternal(_dTime);
								break;
	}
}

void CUnit::InitializeSimplePressure(double _dTime) {
	//Add plot for droplet size distribution
	AddPlot("Probability density function", "Droplet size [mm]", "PDF", "Time");
	AddPlot("Cumulative distribution function", "Droplet size [mm]", "CDF", "Time");
	
	//Add results: calculated Sauter, median and mean diameter
	AddStateVariable("Sauter diameter [mm]", 0, true);
	AddStateVariable("Count median diameter [mm]", 0, true);
	AddStateVariable("Count mean diameter [mm]", 0, true);
	AddStateVariable("Count mode diameter [mm]", 0, true);
	AddStateVariable("Ohnesorge number [-]", 0, true);
}

void CUnit::InitializePneumaticExternal(double _dTime) {
	//Add plot for droplet size distribution
	AddPlot("Probability density function", "Droplet size [mm]", "PDF", "Time");
	AddPlot("Cumulative distribution function", "Droplet size [mm]", "CDF", "Time");

	//Add results: calculated Sauter, median and mean diameter
	AddStateVariable("Sauter diameter [mm]", 0, true);
	AddStateVariable("Count median diameter [mm]", 0, true);
	AddStateVariable("Count mean diameter [mm]", 0, true);
	AddStateVariable("Count mode diameter [mm]", 0, true);
	AddStateVariable("Ohnesorge number [-]", 0, true);
}

void CUnit::InitializePneumaticInternal(double _dTime) {
	//Add plot for droplet size distribution
	AddPlot("Probability density function", "Droplet size [mm]", "PDF", "Time");
	AddPlot("Cumulative distribution function", "Droplet size [mm]", "CDF", "Time");

	//Add results: calculated Sauter, median and mean diameter
	AddStateVariable("Sauter diameter [mm]", 0, true);
	AddStateVariable("Count median diameter [mm]", 0, true);
	AddStateVariable("Count mean diameter [mm]", 0, true);
	AddStateVariable("Count mode diameter [mm]", 0, true);
	AddStateVariable("Ohnesorge number [-]", 0, true);
}

void CUnit::Simulate(double _dTime)
{
	//Simulation for different models
	switch (m_model) {
		case simplePressure:	SimulateSimplePressure(_dTime);
								break;
		case pneumaticExternal:	SimulatePneumaticExternal(_dTime);
								break;
		case pneumaticInternal:	SimulatePneumaticInternal(_dTime);
								break;
	}

	//Save added state variables for the output
	SaveStateVariables(_dTime);
}

void CUnit::SimulateSimplePressure(double _dTime) {

	CMaterialStream* pInStream = GetPortStream("InPort");
	CMaterialStream* pOutStream = GetPortStream("OutPort");

	pOutStream->CopyFromStream(pInStream, _dTime);

	//Define compound IDs: currently this model is only applicable for water/air system
	std::string liquidKey = GetCompoundKey("H2O");
	std::string gasKey = GetCompoundKey("Air");

	/*
	//Check compounds: for further extension with other liquid/gas-systems than water/air
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
	//Sequence of compounds for calculation: 1st is liquid and 2nd is gas. Rearrange the position of liquid and gas compound if needed
	else if (soaComp0 == 2 && soaComp1 == 1) {
		std::string exchange = compounds[0];
		compounds[0] = compounds[1];
		compounds[1] = exchange;
	}
	*/

	//Physical properties of liquid
	double massLiquid = pInStream->GetMassFlow(_dTime);
	double densityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, DENSITY);
	double volumeLiquid = massLiquid / densityLiquid;
	double sigmaLiquid = GetCompoundConstant(liquidKey, CONST_PROP_USER_DEFINED_01);
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, VISCOSITY);

	//Physical properties of gas
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, gasKey, DENSITY);

	//Unit parameters
	double D = GetConstParameterValue("D") * 1e-3; // Convert unit [mm] to [m]
	double deltaP = GetConstParameterValue("liquidDeltaP") * 1e5; // Convert unit [bar] to [Pa]
	double gDev = GetConstParameterValue("geoStdDev");

	//Calculate Sauter-diameter of droplets in [m]
	double sauter = 2.3 * D * pow((4 * volumeLiquid) / (pow(D, 2) * MATH_PI * sqrt(2 * deltaP / densityLiquid)), 0.25) *
							  pow(deltaP * D / sigmaLiquid, -0.25) *
							  pow(D * sqrt(deltaP * densityLiquid) / viscosityLiquid, -0.25) *
							  pow(densityLiquid / densityGas, 0.25);
	//Calculate count median diameter in [m]
	double median = sauter * exp(-2.5 * pow(log(gDev),2));
	//Calculate count mean diameter in [m]
	double mean = median * exp(0.5 * pow(log(gDev), 2));
	//Calculate count mode diameter in [m]
	double mode = median * exp(-0.5 * pow(log(gDev), 2));
	//Calculate Ohnesorge number
	double oh = viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D);

	//Set calculated values for unit output
	SetStateVariable("Sauter diameter [mm]", sauter * 1e3);
	SetStateVariable("Count median diameter [mm]", median * 1e3);
	SetStateVariable("Count mean diameter [mm]", mean * 1e3);
	SetStateVariable("Count mode diameter [mm]", mode * 1e3);
	SetStateVariable("Ohnesorge number [-]", oh);
	
	//Plot PDF and CDF in unit output
	std::vector<double> diameter = GetClassesMeans(DISTR_SIZE);
	std::vector<double> pdf;
	std::vector<double> cdf;
	for (int i = 0; i < diameter.size(); i++) {
		pdf.push_back((1 / (sqrt(2 * MATH_PI)*log(gDev)*diameter[i])) * exp(-0.5 * pow(log(diameter[i] / median) / log(gDev), 2)));
		cdf.push_back(0.5 + 0.5 * erf((log(diameter[i] / median)) / (sqrt(2)*log(gDev))));
	}
	AddCurveOnPlot("Probability density function", _dTime);
	AddPointOnCurve("Probability density function", _dTime, diameter, pdf);
	AddCurveOnPlot("Cumulative distribution function", _dTime);
	AddPointOnCurve("Cumulative distribution function", _dTime, diameter, cdf);

	//Set flow information, PDF and CDF for stream output
	pOutStream->SetMassFlow(_dTime, 0);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_SOLID, massLiquid);
	pOutStream->SetPSD(_dTime, PSD_q0, pdf);
	pOutStream->SetPSD(_dTime, PSD_Q0, cdf);
}

void CUnit::SimulatePneumaticExternal(double _dTime) {
	CMaterialStream* pInStream = GetPortStream("InPort");
	CMaterialStream* pOutStream = GetPortStream("OutPort");

	pOutStream->CopyFromStream(pInStream, _dTime);

	//Define compound IDs: currently this model is only applicable for water/air system
	std::string liquidKey = GetCompoundKey("H2O");
	std::string gasKey = GetCompoundKey("Air");

	//Unit parameters
	double D = GetConstParameterValue("liquidD") * 1e-3; // Convert unit [mm] to [m]
	double deltaP = GetConstParameterValue("externalGasDeltaP") * 1e5; // Convert unit [bar] to [Pa]
	double gDev = GetConstParameterValue("geoStdDev");
	double alr = GetConstParameterValue("ALR");

	//Physical properties of inlet liquid
	double massLiquid = pInStream->GetPhaseMassFlow(_dTime, SOA_SOLID);
	double densityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, DENSITY);
	double volumeLiquid = massLiquid / densityLiquid;
	double sigmaLiquid = GetCompoundConstant(liquidKey, CONST_PROP_USER_DEFINED_01);
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, VISCOSITY);

	//Physical properties of inlet gas
	double massGas = alr * massLiquid;
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, gasKey, DENSITY);

	//Calculate Sauter-diameter of droplets in [mm]
	double sauter = 0.35 * D * pow(deltaP * D / (sigmaLiquid * pow(1 + massGas / massLiquid, 2)), -0.4) *
							   (1 + 2.5 * viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D));
	//Calculate count median diameter in [mm]
	double median = sauter * exp(-2.5 * pow(log(gDev), 2));
	//Calculate mean diameter in [mm]
	double mean = median * exp(0.5 * pow(log(gDev), 2));
	//Calculate count mode diameter in [mm]
	double mode = median * exp(-0.5 * pow(log(gDev), 2));
	//Calculate Ohnesorge number
	double oh = viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D);

	//Set calculated values for the output
	SetStateVariable("Sauter diameter [mm]", sauter * 1e3);
	SetStateVariable("Count median diameter [mm]", median * 1e3);
	SetStateVariable("Count mean diameter [mm]", mean * 1e3);
	SetStateVariable("Count mode diameter [mm]", mode * 1e3);
	SetStateVariable("Ohnesorge number [-]", oh);

	//Plot PDF and CDF in unit output
	std::vector<double> diameter = GetClassesMeans(DISTR_SIZE);
	std::vector<double> pdf;
	std::vector<double> cdf;
	for (int i = 0; i < diameter.size(); i++) {
		pdf.push_back((1 / (sqrt(2 * MATH_PI)*log(gDev)*diameter[i])) * exp(-0.5 * pow(log(diameter[i] / median) / log(gDev), 2)));
		cdf.push_back(0.5 + 0.5 * erf((log(diameter[i] / median)) / (sqrt(2)*log(gDev))));
	}
	AddCurveOnPlot("Probability density function", _dTime);
	AddPointOnCurve("Probability density function", _dTime, diameter, pdf);
	AddCurveOnPlot("Cumulative distribution function", _dTime);
	AddPointOnCurve("Cumulative distribution function", _dTime, diameter, cdf);

	//Set flow information, PDF and CDF for stream output
	pOutStream->SetMassFlow(_dTime, 0);
	pOutStream->SetMassFlow(_dTime, massLiquid + massGas);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_SOLID, massLiquid);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_VAPOR, massGas);
	pOutStream->SetPSD(_dTime, PSD_q0, pdf);
	pOutStream->SetPSD(_dTime, PSD_Q0, cdf);
}

void CUnit::SimulatePneumaticInternal(double _dTime) {
	CMaterialStream* pInStream = GetPortStream("InPort");
	CMaterialStream* pOutStream = GetPortStream("OutPort");

	pOutStream->CopyFromStream(pInStream, _dTime);

	//Define compound IDs: currently this model is only applicable for water/air system
	std::string liquidKey = GetCompoundKey("H2O");
	std::string gasKey = GetCompoundKey("Air");

	//Unit parameters
	double D = GetConstParameterValue("liquidD") * 1e-3; // Convert unit [mm] to [m]
	double deltaP = GetConstParameterValue("internalGasDeltaP") * 1e5; // Convert unit [bar] to [Pa]
	double gDev = GetConstParameterValue("geoStdDev");
	double alr = GetConstParameterValue("ALR");

	//Physical properties of inlet liquid
	double massLiquid = pInStream->GetPhaseMassFlow(_dTime, SOA_SOLID);
	double densityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, DENSITY);
	double volumeLiquid = massLiquid / densityLiquid;
	double sigmaLiquid = GetCompoundConstant(liquidKey, CONST_PROP_USER_DEFINED_01);
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, VISCOSITY);

	//Physical properties of inlet gas
	double massGas = alr * massLiquid;
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, gasKey, DENSITY);

	//Calculate Sauter-diameter of droplets in [mm]
	double sauter = 0.4 * D * pow(deltaP * D / (sigmaLiquid * pow(1 + massGas / massLiquid, 2)), -0.4) *
							  (1 + 0.4 * viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D));
	//Calculate count median diameter in [mm]
	double median = sauter * exp(-2.5 * pow(log(gDev), 2));
	//Calculate mean diameter in [mm]
	double mean = median * exp(0.5 * pow(log(gDev), 2));
	//Calculate count mode diameter in [mm]
	double mode = median * exp(-0.5 * pow(log(gDev), 2));
	//Calculate Ohnesorge number
	double oh = viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D);

	//Set calculated values for the output
	SetStateVariable("Sauter diameter [mm]", sauter * 1e3);
	SetStateVariable("Count median diameter [mm]", median * 1e3);
	SetStateVariable("Count mean diameter [mm]", mean * 1e3);
	SetStateVariable("Count mode diameter [mm]", mode * 1e3);
	SetStateVariable("Ohnesorge number [-]", oh);

	//Plot PDF and CDF in unit output
	std::vector<double> diameter = GetClassesMeans(DISTR_SIZE);
	std::vector<double> pdf;
	std::vector<double> cdf;
	for (int i = 0; i < diameter.size(); i++) {
		pdf.push_back((1 / (sqrt(2 * MATH_PI)*log(gDev)*diameter[i])) * exp(-0.5 * pow(log(diameter[i] / median) / log(gDev), 2)));
		cdf.push_back(0.5 + 0.5 * erf((log(diameter[i] / median)) / (sqrt(2)*log(gDev))));
	}
	AddCurveOnPlot("Probability density function", _dTime);
	AddPointOnCurve("Probability density function", _dTime, diameter, pdf);
	AddCurveOnPlot("Cumulative distribution function", _dTime);
	AddPointOnCurve("Cumulative distribution function", _dTime, diameter, cdf);

	//Set flow information, PDF and CDF for stream output
	pOutStream->SetMassFlow(_dTime, 0);
	pOutStream->SetMassFlow(_dTime, massLiquid + massGas);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_SOLID, massLiquid);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_VAPOR, massGas);
	pOutStream->SetPSD(_dTime, PSD_q0, pdf);
	pOutStream->SetPSD(_dTime, PSD_Q0, cdf);
}

void CUnit::Finalize()
{

}