// EnergyPlus Headers
#include <WindowASHRAE1588RP.hh>
#include <ConvectionCoefficients.hh>
#include <DataGlobals.hh>
#include <DataHeatBalance.hh>
#include <DataHeatBalFanSys.hh>
#include <DataIPShortCuts.hh>
#include <DataSurfaces.hh>
#include <HeatBalanceManager.hh>
#include <InputProcessor.hh>
#include <WindowManager.hh>

namespace EnergyPlus {

namespace WindowASHRAE1588RP {

using namespace DataGlobals;
using namespace DataHeatBalance;
using namespace DataHeatBalFanSys;
using namespace DataIPShortCuts;
using namespace DataSurfaces;

using ConvectionCoefficients::SetExtConvectionCoeff;
using ConvectionCoefficients::CalcISO15099WindowIntConvCoeff;
using InputProcessor::GetNumObjectsFound;
using InputProcessor::GetObjectItem;
using InputProcessor::VerifyName;
using HeatBalanceManager::SetupSimpleWindowGlazingSystem;
using WindowManager::CalcWindowHeatBalance;

Fstring CurrentModuleObject( MaxNameLength ); // to assist in getting input

void
CreateASHRAE1588RPConstructions( int & ConstrNum, bool & ErrorsFound )
{

	// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
	int ConstructNumAlpha; // Number of construction alpha names being passed
	int DummyNumProp; // dummy variable for properties being passed
	int IOStat; // IO Status when calling get input subroutine
	FArray1D_Fstring ConstructAlphas( 1, sFstring( MaxNameLength ) ); // Construction Alpha names defined
	FArray1D< Real64 > ConstructionNums( 3 ); // Temporary array to transfer construction properties
	bool ErrorInName;
	bool IsBlank;
	int Loop;

	int TotWinASHRAE1588Constructs = GetNumObjectsFound( "Construction:WindowASHRAE1588RP" ); // Number of window constructions based on ASHRAE 1588RP

	CurrentModuleObject = "Construction:WindowASHRAE1588RP";
	for ( Loop = 1; Loop <= TotWinASHRAE1588Constructs; ++Loop ) { // Loop through all WindowASHRAE1588RP constructions.

		//Get the object names for each construction from the input processor
		GetObjectItem( CurrentModuleObject, Loop, ConstructAlphas, ConstructNumAlpha, ConstructionNums, DummyNumProp, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

		ErrorInName = false;
		IsBlank = false;
		VerifyName( ConstructAlphas( 1 ), Construct.Name(), ConstrNum, ErrorInName, IsBlank, trim( CurrentModuleObject ) + " Name" );
		if ( IsBlank ) {
			ErrorsFound = true;
			continue;
		}

		++ConstrNum;

		Construct( ConstrNum ).Name = ConstructAlphas( 1 );
		Construct( ConstrNum ).TypeIsWindow = true;

		// set initial guesses

		int number_of_panes = 1;
		int number_of_gaps = number_of_panes - 1;

		// window type variables (initial guess = horizontal slider)
		Real64 width = 1.5;
		Real64 height = 1.2;
		Real64 cos_tilt = 0.0; // 90 deg

		int previous_number_materials = TotMaterials;

		// Allocate temporary arrays
		Surface.allocate(1);
		SurfaceWindow.allocate(1);
		Zone.allocate(1);
		TempEffBulkAir.allocate(1);
		HConvIn.allocate(1);
		MAT.allocate(1);
		QHTRadSysSurf.allocate(1);
		QHWBaseboardSurf.allocate(1);
		QSteamBaseboardSurf.allocate(1);
		QElecBaseboardSurf.allocate(1);

		// Create dummy zone
		int ZoneNum = 1;
		//Zone( ZoneNum )

		// Create dummy Surface to perform heat transfer calcs
		int SurfNum = 1;
		Surface( SurfNum ).Name = ConstructAlphas( 1 ) + ":Surface";
		Surface( SurfNum ).Class = SurfaceClass_Window;
		Surface( SurfNum ).Tilt = 90.0;
		Surface( SurfNum ).CosTilt = cos_tilt;
		Surface( SurfNum ).HeatTransSurf = true;
		// Skip base surface stuff?
		Surface( SurfNum ).ExtBoundCond = 0;
		Surface( SurfNum ).ExtSolar = true;
		Surface( SurfNum ).ExtWind = true;
		Surface( SurfNum ).Zone = 1;
		Surface( SurfNum ).ViewFactorGround = 0.5;
		Surface( SurfNum ).TAirRef = ZoneMeanAirTemp;

		SurfaceWindow( SurfNum ).ShadingFlag = -1;



		// This is where the iterative optimization loop will begin
		int number_of_new_materials = number_of_panes + number_of_gaps;

		TotMaterials += number_of_new_materials;

		// Create New Material objects
		MaterialSave.allocate( previous_number_materials );
		NominalRSave.allocate( previous_number_materials );
		MaterialSave( {1,previous_number_materials} ) = Material( {1,previous_number_materials} );
		NominalRSave( {1,previous_number_materials} ) = NominalR( {1,previous_number_materials} );
		Material.deallocate();
		NominalR.deallocate();
		Material.allocate( TotMaterials );
		NominalR.allocate( TotMaterials );
		Material( {1,TotMaterials - number_of_new_materials} ) = MaterialSave( {1,TotMaterials - number_of_new_materials} );
		NominalR( {1,TotMaterials - number_of_new_materials} ) = NominalRSave( {1,TotMaterials - number_of_new_materials} );
		MaterialSave.deallocate();
		NominalRSave.deallocate();

		// Define material properties (currently this is the same as the simple glazing system properties)
		Material( TotMaterials ).Group = WindowSimpleGlazing;
		Material( TotMaterials ).Name = ConstructAlphas( 1 ) + ":Pane1";
		Material( TotMaterials ).SimpleWindowUfactor = ConstructionNums( 1 );
		Material( TotMaterials ).SimpleWindowSHGC = ConstructionNums( 2 );
		if ( ! lNumericFieldBlanks( 3 ) ) {
			Material( TotMaterials ).SimpleWindowVisTran = ConstructionNums( 3 );
			Material( TotMaterials ).SimpleWindowVTinputByUser = true;
		}


		SetupSimpleWindowGlazingSystem( TotMaterials );

		Construct( ConstrNum ).TotLayers = 1;
		Construct( ConstrNum ).LayerPoint( 1 ) = TotMaterials;
		Construct( ConstrNum ).TotGlassLayers = number_of_panes;

		for ( int Layer = 1; Layer <= Construct( ConstrNum ).TotLayers; ++Layer ) {
			NominalRforNominalUCalculation( ConstrNum ) += NominalR( Construct( ConstrNum ).LayerPoint( Layer ) );
		}

		Surface( SurfNum ).Construction = ConstrNum;

		// TODO Create new WindowFrameAndDivider objects (see Window 5 method)

		// Calculate U-factor and SHGC and compare to inputs

		// U-factor conditions
		Real64 air_temp = 21.0;

		Surface( SurfNum ).OutDryBulbTemp = -18.0;
		SurfaceWindow( SurfNum ).IRfromParentZone = StefanBoltzmann*std::pow(air_temp + KelvinConv,4);

		Real64 inside_temp = air_temp;
		Real64 outside_temp = Surface( SurfNum ).OutDryBulbTemp; // initial guess
//		Real64 h_exterior = SetExtConvectionCoeff( SurfNum );
//		HConvIn( SurfNum ) = CalcISO15099WindowIntConvCoeff( SurfNum, inside_temp, air_temp);
//		CalcWindowHeatBalance( SurfNum, h_exterior, inside_temp, outside_temp );

		// if match not obtained adjust inputs

		// end loop

		// deallocate temporary arrays
		Surface.deallocate();
		SurfaceWindow.deallocate();
		Zone.deallocate();
		TempEffBulkAir.deallocate();
		HConvIn.deallocate();
		MAT.deallocate();
		QHTRadSysSurf.deallocate();
		QHWBaseboardSurf.deallocate();
		QSteamBaseboardSurf.deallocate();
		QElecBaseboardSurf.deallocate();


	} // ...end of WindowASHRAE1588RP Constructions DO loop


}


} // WindowASHRAE1588RP

} // EnergyPlus
