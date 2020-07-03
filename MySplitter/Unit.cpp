#define DLL_EXPORT
#include "Unit.h"

extern "C" DECLDIR CBaseUnit* DYSSOL_CREATE_MODEL_FUN()
{
	return new CUnit();
}

//////////////////////////////////////////////////////////////////////////
/// A splitter with one inlet and three outlets

CUnit::CUnit()
{
	/// Basic unit's info ///
	m_sUnitName = "mySplitter";
	m_sAuthorName = "XYZhou";
	m_sUniqueID = "80E3E38A-98E4-4D1A-AD07-8D20A54B8B58";

	/// Add ports ///
	AddPort("In", INPUT_PORT);
	AddPort("Out1", OUTPUT_PORT);
	AddPort("Out2", OUTPUT_PORT);
	AddPort("Out3", OUTPUT_PORT);

	/// Add unit parameters ///
	AddConstParameter("k1", 0, 1, 0, "k1", "splitting factor for outlet1");
	AddConstParameter("k2", 0, 1, 0, "k2", "splitting factor for outlet2");
}

CUnit::~CUnit()
{

}

void CUnit::Initialize(double _dTime)
{
	/// Add state variables ///
	/// AddStateVariable("VarName", 0, true);

}

void CUnit::Simulate(double _dTime)
{
	CMaterialStream* pInStream = GetPortStream("In");
	CMaterialStream* pOutStream1 = GetPortStream("Out1");
	CMaterialStream* pOutStream2 = GetPortStream("Out2");
	CMaterialStream* pOutStream3 = GetPortStream("Out3");

	/// Copy input stream information to output streams
	pOutStream1->CopyFromStream(pInStream, _dTime);
	pOutStream2->CopyFromStream(pInStream, _dTime);
	pOutStream3->CopyFromStream(pInStream, _dTime);

	/// Calculate the outlets
	double InletMassFlow = pInStream->GetMassFlow(_dTime);
	double SplittFactor1 = GetConstParameterValue("k1");
	double SplittFactor2 = GetConstParameterValue("k2");
	if (SplittFactor1 + SplittFactor2 > 1) {
		RaiseError("3rd outlet is minus! Please ensure the sum of both splitt factors is not greater than 1.");
	}
	double outletMassFlow1 = InletMassFlow * SplittFactor1;
	double outletMassFlow2 = InletMassFlow * SplittFactor2;
	double outletMassFlow3 = InletMassFlow * (1 - SplittFactor1 - SplittFactor2);
	pOutStream1->SetMassFlow(_dTime, outletMassFlow1, BASIS_MASS);
	pOutStream2->SetMassFlow(_dTime, outletMassFlow2, BASIS_MASS);
	pOutStream3->SetMassFlow(_dTime, outletMassFlow3, BASIS_MASS);
}

void CUnit::Finalize()
{

}