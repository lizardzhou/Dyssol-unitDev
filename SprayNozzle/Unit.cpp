#define DLL_EXPORT
#include "Unit.h"
#include "math.h"
#include <algorithm>

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// Spray nozzle with 4 different models: 
///  1) single-fluid; 
///  2) two-fluid with external mixing; 
///  3) two-fluid with internal mixing; 
///  4) Rotary;
/// User only needs to define the liquid phase as solid. In case of pneumatic model, the gas phase will be calculated using ALR
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
	AddGroupParameter("Model", singleFluid, { singleFluid, twoFluidExternal, twoFluidInternal, rotary }, { "Simple pressure nozzle", "Two-fluid nozzle with external mixing", "Two-fluid nozzle with internal mixing", "Rotary nozzle" }, "Atomization model");
	AddConstParameter("D", 1, 2, 1, "mm", "Liquid oriface diameter"); // nozzle diameter for single-fliud model
	AddConstParameter("liquidD", 0.5, 2, 1, "mm", "Liquid oriface diameter"); // nozzle diameter for two-fluid model
	AddConstParameter("diskD", 40, 120, 50, "mm", "Disk diameter"); // disk diameter for rotary model
	AddConstParameter("R", 40, 180, 50, "mm", "Downstream distance along the spray trajectory"); // downstream distance along the spray trajectory
	AddConstParameter("omega", 830, 1700, 1500, "rad/s", "Angular velocity"); // angular speed for rotary model
	AddConstParameter("ALR", 1e-5, 5, 1, "-", "Air-to-liquid mass flow ratio"); // air-to-liquid mass flow ratio for two-fluid nozzle
	AddConstParameter("liquidDeltaP", 3, 200, 50, "bar", "Pressure drop in nozzle"); //liquid pressure drop in single-fluid nozzle
	AddConstParameter("externalGasDeltaP", 0.5, 3, 3, "bar", "Gas pressure drop in nozzle"); // gas pressure drop in two-fluid nozzle with external mixing
	AddConstParameter("internalGasDeltaP", 0.5, 15, 3, "bar", "Gas pressure drop in nozzle"); // gas pressure drop in two-fluid nozzle with internal mixing
	AddConstParameter("geoStdDev", 1, 5, 1.5, "-", "Geometric standard deviation of droplet size distribution"); // geometric standard deviation
	AddStringParameter("Drop compound", "H2O", "Show the drop compound name");
	AddParametersToGroup("Model", "Simple pressure nozzle", { "D", "liquidDeltaP", "geoStdDev" });
	AddParametersToGroup("Model", "Two-fluid nozzle with external mixing", { "liquidD", "externalGasDeltaP", "ALR", "geoStdDev" });
	AddParametersToGroup("Model", "Two-fluid nozzle with internal mixing", { "liquidD", "internalGasDeltaP", "ALR", "geoStdDev" });
	AddParametersToGroup("Model", "Rotary nozzle", {"diskD", "R", "omega", "geoStdDev"});
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

void CUnit::calcDiameter(double &sauter, double &gDev, double &median, double &mean, double &mode) {
	//Calculate count median diameter in [mm]
	median = sauter * exp(-2.5 * pow(log(gDev), 2));
	//Calculate mean diameter in [mm]
	mean = median * exp(0.5 * pow(log(gDev), 2));
	//Calculate count mode diameter in [mm]
	mode = median * exp(-pow(log(gDev), 2));
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

	/*
	//Initialization for different models
	switch (m_model) {
	case singleFluid:	InitializeOhnesorge(_dTime);
	break;
	case twoFluidExternal:	InitializeOhnesorge(_dTime);
	break;
	case twoFluidInternal: InitializeOhnesorge(_dTime);
	break;
	case rotary:			InitializeWeber(_dTime);
	break;
	*/

	//Add plot for droplet size distribution
	AddPlot("Probability density function", "Droplet size [mm]", "PDF", "Time");

	//Add results: calculated Sauter, median and mean diameter, dimensionless numbers
	AddStateVariable("Sauter diameter [micron]", 0, true);
	AddStateVariable("Count median diameter [micron]", 0, true);
	AddStateVariable("Count mean diameter [micron]", 0, true);
	AddStateVariable("Count mode diameter [micron]", 0, true);
	AddStateVariable("Droplet Ohnesorge number [-]", 0, true);
	AddStateVariable("Droplet Weber number [-]", 0, true);
}

/*
void CUnit::InitializeOhnesorge(double _dTime) {
	//Add plot for droplet size distribution
	AddPlot("Probability density function", "Droplet size [mm]", "PDF", "Time");
	//AddPlot("Cumulative distribution function", "Droplet size [mm]", "CDF", "Time");
	
	//Add results: calculated Sauter, median and mean diameter, Ohnesorge number
	AddStateVariable("Sauter diameter [mm]", 0, true);
	AddStateVariable("Count median diameter [mm]", 0, true);
	AddStateVariable("Count mean diameter [mm]", 0, true);
	AddStateVariable("Count mode diameter [mm]", 0, true);
	AddStateVariable("Ohnesorge number [-]", 0, true);
	AddStateVariable("Weber number [-]", 0, true);
}
*/

/*
void CUnit::InitializeWeber(double _dTime) {
	//Add plot for droplet size distribution
	AddPlot("Probability density function", "Droplet size [mm]", "PDF", "Time");
	//AddPlot("Cumulative distribution function", "Droplet size [mm]", "CDF", "Time");

	//Add results: calculated Sauter, median and mean diameter, Weber number
	AddStateVariable("Sauter diameter [mm]", 0, true);
	AddStateVariable("Count median diameter [mm]", 0, true);
	AddStateVariable("Count mean diameter [mm]", 0, true);
	AddStateVariable("Count mode diameter [mm]", 0, true);
	AddStateVariable("Weber number [-]", 0, true);
}
*/

void CUnit::Simulate(double _dTime)
{
	//Simulation for different models
	switch (m_model) {
		case singleFluid:	SimulateSingleFluid(_dTime);
		break;
		case twoFluidExternal:	SimulateTwoFluidExternal(_dTime);
		break;
		case twoFluidInternal:	SimulateTwoFluidInternal(_dTime);
		break;
		case rotary:			SimulateRotary(_dTime);
		break;
	}

	//Save added state variables for the output
	SaveStateVariables(_dTime);
}

void CUnit::SimulateSingleFluid(double _dTime) {

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
	//Calculate count median, mean and mode diameter in [m]
	double median = 0;
	double mean = 0;
	double mode = 0;
	calcDiameter(sauter, gDev, median, mean, mode);
	//Calculate Ohnesorge number
	double ohD = viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * sauter);
	//Calculate Weber number
	double velocityLiquid = volumeLiquid / (MATH_PI * pow(D, 2) / 4);
	double weD = densityGas * pow(velocityLiquid, 2) * sauter / sigmaLiquid;

	//Set calculated values for unit output
	SetStateVariable("Sauter diameter [micron]", sauter * 1e6);
	SetStateVariable("Count median diameter [micron]", median * 1e6);
	SetStateVariable("Count mean diameter [micron]", mean * 1e6);
	SetStateVariable("Count mode diameter [micron]", mode * 1e6);
	SetStateVariable("Droplet Ohnesorge number [-]", ohD);
	SetStateVariable("Droplet Weber number [-]", weD);
	
	//Plot PDF and CDF in unit output
	std::vector<double> diameter = GetClassesMeans(DISTR_SIZE);
	std::vector<double> pdf;
	for (int i = 0; i < diameter.size(); i++) {
		pdf.push_back((1 / (sqrt(2 * MATH_PI)*log(gDev)*diameter[i])) * exp(-0.5 * pow(log(diameter[i] / median) / log(gDev), 2)));
		//extra calculation of CFD is not needed because Dyssol will do this internally
	}
	AddCurveOnPlot("Probability density function", _dTime);
	AddPointOnCurve("Probability density function", _dTime, diameter, pdf);

	//Check if a PDF-plot is incomplete
	double lastPDF = *(pdf.end() - 1); // last element in PDF
	double maxPDF = *std::max_element(pdf.begin(), pdf.end()); // largest element in PDF
	if (!(lastPDF < 1e-7 && maxPDF > 1)) { 
		RaiseWarning("The defined input distribution range is too short to display the complete PDF, leading to wrong results!");
	}
	//maxDRange < 1e-7: the PDF-value at the end should approch 0; 
	//maxD > 1: if the defined input range is so narrow that all PDF-values are starting values approaching 0, in this case the curve is also incomplete

	//Set flow information, PDF and CDF for stream output
	pOutStream->SetMassFlow(_dTime, 0);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_SOLID, massLiquid);
	pOutStream->SetPSD(_dTime, PSD_q0, pdf);
}

void CUnit::SimulateTwoFluidExternal(double _dTime) {
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

	//Calculate Sauter-diameter of droplets in [m]
	double sauter = 0.35 * D * pow(deltaP * D / (sigmaLiquid * pow(1 + massLiquid / massGas, 2)), -0.4) *
							   (1 + 2.5 * viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D));
	
	//Calculate count median, mean and mode diameter in [m]
	double median = 0;
	double mean = 0;
	double mode = 0;
	calcDiameter(sauter, gDev, median, mean, mode);
	//Calculate Ohnesorge number
	double ohD = viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * sauter);
	//Calculate Weber number
	double velocityLiquid = volumeLiquid / (MATH_PI * pow(D, 2) / 4);
	double weD = densityGas * pow(velocityLiquid, 2) * sauter / sigmaLiquid;

	//Set calculated values for the output
	SetStateVariable("Sauter diameter [micron]", sauter * 1e6);
	SetStateVariable("Count median diameter [micron]", median * 1e6);
	SetStateVariable("Count mean diameter [micron]", mean * 1e6);
	SetStateVariable("Count mode diameter [micron]", mode * 1e6);
	SetStateVariable("Droplet Ohnesorge number [-]", ohD);
	SetStateVariable("Droplet Weber number [-]", weD);

	//Plot PDF and CDF in unit output
	std::vector<double> diameter = GetClassesMeans(DISTR_SIZE);
	std::vector<double> pdf;
	//std::vector<double> cdf;
	for (int i = 0; i < diameter.size(); i++) {
		pdf.push_back((1 / (sqrt(2 * MATH_PI)*log(gDev)*diameter[i])) * exp(-0.5 * pow(log(diameter[i] / median) / log(gDev), 2)));
	}
	AddCurveOnPlot("Probability density function", _dTime);
	AddPointOnCurve("Probability density function", _dTime, diameter, pdf);

	//Check if a PDF-plot is incomplete
	double lastPDF = *(pdf.end() - 1); // last element in PDF
	double maxPDF = *std::max_element(pdf.begin(), pdf.end()); // largest element in PDF
	if (!(lastPDF < 1e-7 && maxPDF > 1)) {
		RaiseWarning("The defined input distribution range is too short to display the complete PDF, leading to wrong results!");
	}
	//maxDRange < 1e-7: the PDF-value at the end should approch 0; 
	//maxD > 1: if the defined input range is so narrow that all PDF-values are starting values approaching 0, in this case the curve is also incomplete

	//Set flow information, PDF and CDF for stream output
	pOutStream->SetMassFlow(_dTime, 0);
	pOutStream->SetMassFlow(_dTime, massLiquid + massGas);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_SOLID, massLiquid);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_VAPOR, massGas);
	pOutStream->SetPSD(_dTime, PSD_q0, pdf);
}

void CUnit::SimulateTwoFluidInternal(double _dTime) {
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

	//Calculate Sauter-diameter of droplets in [m]
	double sauter = 0.4 * D * pow(deltaP * D / (sigmaLiquid * pow(1 + massLiquid / massGas, 2)), -0.4) *
							  (1 + 0.4 * viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * D));
	//Calculate count median, mean and mode diameter in [m]
	double median = 0;
	double mean = 0;
	double mode = 0;
	calcDiameter(sauter, gDev, median, mean, mode);
	//Calculate Ohnesorge numbers
	double ohD = viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * sauter);
	//Calculate Weber number
	double velocityLiquid = volumeLiquid / (MATH_PI * pow(D, 2) / 4);
	double weD = densityGas * pow(velocityLiquid, 2) * sauter / sigmaLiquid;

	//Set calculated values for the output
	SetStateVariable("Sauter diameter [micron]", sauter * 1e6);
	SetStateVariable("Count median diameter [micron]", median * 1e6);
	SetStateVariable("Count mean diameter [micron]", mean * 1e6);
	SetStateVariable("Count mode diameter [micron]", mode * 1e6);
	SetStateVariable("Droplet Ohnesorge number [-]", ohD);
	SetStateVariable("Droplet Weber number [-]", weD);

	//Plot PDF and CDF in unit output
	std::vector<double> diameter = GetClassesMeans(DISTR_SIZE);
	std::vector<double> pdf;
	for (int i = 0; i < diameter.size(); i++) {
		pdf.push_back((1 / (sqrt(2 * MATH_PI)*log(gDev)*diameter[i])) * exp(-0.5 * pow(log(diameter[i] / median) / log(gDev), 2)));
	}
	AddCurveOnPlot("Probability density function", _dTime);
	AddPointOnCurve("Probability density function", _dTime, diameter, pdf);

	//Check if a PDF-plot is incomplete
	double lastPDF = *(pdf.end() - 1); // last element in PDF
	double maxPDF = *std::max_element(pdf.begin(), pdf.end()); // largest element in PDF
	if (!(lastPDF < 1e-7 && maxPDF > 1)) {
		RaiseWarning("The defined input distribution range is too short to display the complete PDF, leading to wrong results!");
	}
	//maxDRange < 1e-7: the PDF-value at the end should approch 0; 
	//maxD > 1: if the defined input range is so narrow that all PDF-values are starting values approaching 0, in this case the curve is also incomplete

	//Set flow information, PDF and CDF for stream output
	pOutStream->SetMassFlow(_dTime, 0);
	pOutStream->SetMassFlow(_dTime, massLiquid + massGas);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_SOLID, massLiquid);
	pOutStream->SetPhaseMassFlow(_dTime, SOA_VAPOR, massGas);
	pOutStream->SetPSD(_dTime, PSD_q0, pdf);
}

void CUnit::SimulateRotary(double _dTime) {
	CMaterialStream* pInStream = GetPortStream("InPort");
	CMaterialStream* pOutStream = GetPortStream("OutPort");

	pOutStream->CopyFromStream(pInStream, _dTime);

	//Define compound IDs: currently this model is only applicable for water/air system
	std::string liquidKey = GetCompoundKey("H2O");
	std::string gasKey = GetCompoundKey("Air");

	//Unit parameters
	double D = GetConstParameterValue("diskD") * 1e-3; // Convert unit [mm] to [m]
	double R = GetConstParameterValue("R") * 1e-3; // Convert unit [mm] to [m]
	double omega = GetConstParameterValue("omega") / (2 * MATH_PI);
	double gDev = GetConstParameterValue("geoStdDev");

	//Physical properties of inlet liquid
	double massLiquid = pInStream->GetPhaseMassFlow(_dTime, SOA_SOLID);
	double densityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, DENSITY);
	double volumeLiquid = massLiquid / densityLiquid;
	double sigmaLiquid = GetCompoundConstant(liquidKey, CONST_PROP_USER_DEFINED_01);
	double viscosityLiquid = pInStream->GetCompoundTPDProp(_dTime, liquidKey, VISCOSITY);
	
	//Physical properties of inlet gas
	double densityGas = pInStream->GetCompoundTPDProp(_dTime, gasKey, DENSITY);

	//Calculate Sauter-diameter of droplets in [m]
	double sauter = 27.81 * D * pow(volumeLiquid / (omega*pow(D, 3)), 0.051) * pow(R/D, 0.581) *
								pow(densityLiquid * omega * pow(D, 2) / viscosityLiquid, -0.651) * 
								pow(pow(D, 3) * pow(omega, 2) * densityLiquid / sigmaLiquid, -0.0218);
	//Calculate count median, mean and mode diameter in [m]
	double median = 0;
	double mean = 0;
	double mode = 0;
	calcDiameter(sauter, gDev, median, mean, mode);
	//Calculate Ohnesorge number
	double ohD = viscosityLiquid / sqrt(sigmaLiquid * densityLiquid * sauter);
	//Calculate Weber numbers
	double velocityLiquid = omega * D;
	double weD = densityGas * pow(velocityLiquid, 2) * sauter / sigmaLiquid;

	//Set calculated values for the output
	SetStateVariable("Sauter diameter [micron]", sauter * 1e6);
	SetStateVariable("Count median diameter [micron]", median * 1e6);
	SetStateVariable("Count mean diameter [micron]", mean * 1e6);
	SetStateVariable("Count mode diameter [micron]", mode * 1e6);
	SetStateVariable("Droplet Ohnesorge number [-]", ohD);
	SetStateVariable("Droplet Weber number [-]", weD);

	//Plot PDF and CDF in unit output
	std::vector<double> diameter = GetClassesMeans(DISTR_SIZE);
	std::vector<double> pdf;
	for (int i = 0; i < diameter.size(); i++) {
		pdf.push_back((1 / (sqrt(2 * MATH_PI)*log(gDev)*diameter[i])) * exp(-0.5 * pow(log(diameter[i] / median) / log(gDev), 2)));
	}
	AddCurveOnPlot("Probability density function", _dTime);
	AddPointOnCurve("Probability density function", _dTime, diameter, pdf);

	//Check if a PDF-plot is incomplete
	double lastPDF = *(pdf.end() - 1); // last element in PDF
	double maxPDF = *std::max_element(pdf.begin(), pdf.end()); // largest element in PDF
	if (!(lastPDF < 1e-7 && maxPDF > 1)) {
		RaiseWarning("The defined input distribution range is too short to display the complete PDF, leading to wrong results!");
	}
	//maxDRange < 1e-7: the PDF-value at the end should approch 0; 
	//maxD > 1: if the defined input range is so narrow that all PDF-values are starting values approaching 0, in this case the curve is also incomplete

	//Set flow information, PDF and CDF for stream output
	pOutStream->SetMassFlow(_dTime, 0);
	pOutStream->SetMassFlow(_dTime, massLiquid);
	pOutStream->SetPSD(_dTime, PSD_q0, pdf);
}

void CUnit::Finalize()
{

}