// EnergyPlus Headers
#include <WindowASHRAE1588RP.hh>
#include <ConvectionCoefficients.hh>
#include <DataEnvironment.hh>
#include <DataGlobals.hh>
#include <DataHeatBalance.hh>
#include <DataHeatBalFanSys.hh>
#include <DataHeatBalSurface.hh>
#include <DataIPShortCuts.hh>
#include <DataSurfaces.hh>
#include <HeatBalanceManager.hh>
#include <HeatBalanceSurfaceManager.hh>
#include <InputProcessor.hh>
#include <WindowManager.hh>
#include <SolarShading.hh>

namespace EnergyPlus {

namespace WindowASHRAE1588RP {

using namespace DataEnvironment;
using namespace DataGlobals;
using namespace DataHeatBalance;
using namespace DataHeatBalFanSys;
using namespace DataHeatBalSurface;
using namespace DataIPShortCuts;
using namespace DataSurfaces;

using ConvectionCoefficients::SetExtConvectionCoeff;
using ConvectionCoefficients::CalcISO15099WindowIntConvCoeff;
using InputProcessor::GetNumObjectsFound;
using InputProcessor::GetObjectItem;
using InputProcessor::VerifyName;
using HeatBalanceManager::SetupSimpleWindowGlazingSystem;
using HeatBalanceSurfaceManager::InitSolarHeatGains;
using WindowManager::CalcWindowHeatBalance;
using WindowManager::InitGlassOpticalCalculations;
using SolarShading::CalcInteriorSolarDistribution;
using SolarShading::ISABSF;

std::string CurrentModuleObject; // to assist in getting input

void
CreateASHRAE1588RPConstructions( int & ConstrNum, bool & ErrorsFound )
{

	int ConstructNumAlpha; // Number of construction alpha names being passed
	int ConstructNumNumeric; // dummy variable for properties being passed
	int IOStat; // IO Status when calling get input subroutine
	FArray1D_string ConstructAlphas( 1 ); // Construction Alpha names defined
	FArray1D< Real64 > ConstructNumerics( 3 ); // Temporary array to transfer construction properties
	bool ErrorInName;
	bool IsBlank;
	int Loop;

	int TotWinASHRAE1588Constructs = GetNumObjectsFound( "Construction:WindowASHRAE1588RP" ); // Number of window constructions based on ASHRAE 1588RP

	CurrentModuleObject = "Construction:WindowASHRAE1588RP";
	for ( Loop = 1; Loop <= TotWinASHRAE1588Constructs; ++Loop ) { // Loop through all WindowASHRAE1588RP constructions.

		//Get the object names for each construction from the input processor
		GetObjectItem( CurrentModuleObject, Loop, ConstructAlphas, ConstructNumAlpha, ConstructNumerics, ConstructNumNumeric, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

		ErrorInName = false;
		IsBlank = false;
		VerifyName( ConstructAlphas( 1 ), Construct.Name(), ConstrNum, ErrorInName, IsBlank, trim( CurrentModuleObject ) + " Name" );
		if ( IsBlank ) {
			ErrorsFound = true;
			continue;
		}

		++ConstrNum;

		// Save Constructions -- The list will be deleted so that the only
		// construction is the one currently being set for any borrowed subroutines
		FArray1D< Real64 > NominalRforNominalUCalculationSave;

		ConstructSave.allocate( TotConstructs );
		ConstructSave = Construct;
		NominalRforNominalUCalculationSave.allocate( TotConstructs );
		NominalUSave.allocate( TotConstructs );
		NominalRforNominalUCalculationSave = NominalRforNominalUCalculation;
		NominalUSave = NominalU;

		int TotConstructsSave = TotConstructs;

		Construct.deallocate();
		NominalRforNominalUCalculation.deallocate();
		NominalU.deallocate();

		Construct.allocate(1);
		NominalRforNominalUCalculation.allocate(1);
		NominalU.allocate(1);

		TotConstructs = 1;

		ConstructionData new_construct;

		new_construct.Name = ConstructAlphas( 1 );
		new_construct.TypeIsWindow = true;


		// set initial guesses

		int number_of_panes = 1;
		int number_of_gaps = number_of_panes - 1;

		// window type variables (initial guess = horizontal slider)
		Real64 width = 1.5;
		Real64 height = 1.2;
		Real64 tilt = Pi/2; // 90 deg

		int previous_number_materials = TotMaterials;

		// Allocate temporary arrays
		create_dummy_variables();

		Surface( 1 ).Name = ConstructAlphas( 1 ) + ":Surface";
		Surface( 1 ).Tilt = tilt*180/Pi;
		Surface( 1 ).CosTilt = cos(tilt);
		Surface( 1 ).SinTilt = sin(tilt);
		Surface( 1 ).Height = height;
		Surface( 1 ).Area = height*width;
		Surface( 1 ).ViewFactorSky = 0.5 * ( 1.0 + Surface( 1 ).CosTilt );
		Surface( 1 ).ViewFactorGround = 0.5 * ( 1.0 - Surface( 1 ).CosTilt );
		Surface( 1 ).ViewFactorSkyIR = Surface( 1 ).ViewFactorSky;
		Surface( 1 ).ViewFactorGroundIR = Surface( 1 ).ViewFactorGround;
		AirSkyRadSplit( 1 ) = std::sqrt( 0.5 * ( 1.0 + Surface( 1 ).CosTilt ) );

		ASHRAE1588RP_Flag = true;
		KickOffSimulation = false;

		// This is where the iterative optimization loop will begin
		int number_of_new_materials = number_of_panes + number_of_gaps;

		TotMaterials += number_of_new_materials;

		// Construction specific allocations
		AWinSurf.allocate(1, number_of_panes);
		QRadSWwinAbs.allocate(1, number_of_panes);
		QRadSWwinAbsLayer.allocate(1, number_of_panes);

		// Create New Material objects
		// TODO We'll probably need to save the materials and just pass the few materials to the global list that we need. --Similar to constructions
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
		Material( TotMaterials ).SimpleWindowUfactor = ConstructNumerics( 1 );
		Material( TotMaterials ).SimpleWindowSHGC = ConstructNumerics( 2 );
		if ( ! lNumericFieldBlanks( 3 ) ) {
			Material( TotMaterials ).SimpleWindowVisTran = ConstructNumerics( 3 );
			Material( TotMaterials ).SimpleWindowVTinputByUser = true;
		}


		SetupSimpleWindowGlazingSystem( TotMaterials );

		new_construct.TotLayers = 1;
		new_construct.LayerPoint( 1 ) = TotMaterials;
		new_construct.TotGlassLayers = number_of_panes;

		Construct( 1 ) = new_construct;

		for ( int Layer = 1; Layer <= Construct( 1 ).TotLayers; ++Layer ) {
			NominalRforNominalUCalculation( 1 ) += NominalR( Construct( 1 ).LayerPoint( Layer ) );
		}

		Surface( 1 ).Construction = 1; // This is the only construction available to the dummy surface. The actual surface will reference the real construction.

		// TODO Create new WindowFrameAndDivider objects (see Window 5 method)

		// Setup functions

		InitGlassOpticalCalculations();

		// Set up U-factor conditions
		Real64 in_air_temp = 21.0;
		Real64 out_air_temp = -18.0;
		Real64 wind_speed = 5.5;
		Real64 solar_inccident = 0.0;

		calc_window_performance(in_air_temp, out_air_temp, wind_speed, solar_inccident);

		Real64 u_factor = -WinHeatGain(1)/(Surface( 1 ).Area*(in_air_temp - out_air_temp));

		// Set up SHGC conditions
		in_air_temp = 24.0;
		out_air_temp = 32.0;
		wind_speed = 2.75;
		solar_inccident = 783.0;

		calc_window_performance(in_air_temp, out_air_temp, wind_speed, solar_inccident);

		Real64 q_total = WinHeatGain(1);

		// NFRC 201-2014 Equation 8-7
		Real64 q_U = u_factor*Surface( 1 ).Area*(out_air_temp - in_air_temp);

		// NFRC 201-2014 Equation 8-2
		Real64 shgc = (q_total - q_U)/(Surface( 1 ).Area*solar_inccident);

		// if match not obtained adjust inputs

		// Deallocate construction specific arrays
		AWinSurf.deallocate();
		QRadSWwinAbs.deallocate();
		QRadSWwinAbsLayer.deallocate();

		// end loop

		ASHRAE1588RP_Flag = false;
		KickOffSimulation = true;

		// deallocate temporary arrays
		remove_dummy_variables();

		// Restore construction list and copy in new construction
		Real64 newU = NominalU( 1 );
		Real64 newR = NominalRforNominalUCalculation( 1 );

		Construct.deallocate();
		NominalRforNominalUCalculation.deallocate();
		NominalU.deallocate();

		TotConstructs = TotConstructsSave;

		Construct.allocate( TotConstructs );
		Construct = ConstructSave;
		NominalRforNominalUCalculation.allocate( TotConstructs );
		NominalU.allocate( TotConstructs );
		NominalRforNominalUCalculation = NominalRforNominalUCalculationSave;
		NominalUSave = NominalU;

		Construct( ConstrNum ) = new_construct;
		NominalRforNominalUCalculation( ConstrNum ) = newR;
		NominalU( ConstrNum ) = newU;

		ConstructSave.deallocate();
		NominalRforNominalUCalculationSave.deallocate();
		NominalUSave.deallocate();


	} // ...end of WindowASHRAE1588RP Constructions DO loop


}

void calc_window_performance(Real64 T_in, Real64 T_out, Real64 v_ws, Real64 I_s)
{
	// Calculate window performance
	Surface( 1 ).OutDryBulbTemp = T_out;
	TempEffBulkAir( 1 ) = T_in;

	SurfaceWindow( 1 ).IRfromParentZone = StefanBoltzmann*std::pow(T_in + KelvinConv,4);

	// initial guess temperatures
	int num_temps = 2 + 2*Construct(1).TotGlassLayers;
	Real64 in_surf_temp = T_in - (1.0/(num_temps-1))*(T_in - T_out);
	Real64 out_surf_temp = T_out + (1.0/(num_temps-1))*(T_in - T_out);

	Real64 h_exterior_f = 4 + v_ws*4;
	Real64 h_exterior;

	BeamSolarRad = I_s;

	if ( I_s > 0.0 ) {
		SunIsUp = true;
	}

	InitSolarHeatGains();
	CalcInteriorSolarDistribution();

	// Calculate heat balance (iteratively solve for surface temperatures)
	Real64 out_surf_temp_prev = out_surf_temp;
	Real64 in_surf_temp_prev = in_surf_temp;

	Real64 out_surf_temp_diff;
	Real64 in_surf_temp_diff;

	int max_iterations = 20;
	Real64 tolerance = 0.1; // deg C

	for (int i = 0; i < max_iterations; i++) {
		CalcISO15099WindowIntConvCoeff( 1, out_surf_temp, T_out); // This subroutine sets the global HConvIn( 1 ) variable. We will use it to set the exterior natural convection.
		h_exterior = h_exterior_f + HConvIn( 1 ); // add natural convection
		CalcISO15099WindowIntConvCoeff( 1, in_surf_temp, T_in); // This time it's actually being used as intended. HConvIn( 1 ) is referenced from the actual heat balance calculation.
		CalcWindowHeatBalance( 1, h_exterior, in_surf_temp, out_surf_temp );

		out_surf_temp_diff = std::fabs(out_surf_temp - out_surf_temp_prev);
		in_surf_temp_diff = std::fabs(in_surf_temp - in_surf_temp_prev);

		if ( (out_surf_temp_diff < tolerance) && (in_surf_temp_diff < tolerance) ) {
			break;
		}

		out_surf_temp_prev = out_surf_temp;
		in_surf_temp_prev = in_surf_temp;

	}




}

void create_dummy_variables()
{

	// Zone
	Zone.allocate(1);
	NumOfZones = 1;
	Zone( 1 ).SurfaceFirst = 1;
	Zone( 1 ).SurfaceLast = 1;

	MAT.allocate(1);
	ZoneAirHumRatAvg.allocate(1);
	ZoneAirHumRat.allocate(1);
	DSZone.allocate(1);
	DGZone.allocate(1);
	DBZoneSSG.allocate(1);
	DBZone.allocate(1);
	ZoneTransSolar.allocate(1);
	ZoneTransSolarEnergy.allocate(1);
	ZoneBmSolFrExtWinsRep.allocate(1);
	ZoneDifSolFrExtWinsRep.allocate(1);
	ZoneBmSolFrExtWinsRepEnergy.allocate(1);
	ZoneDifSolFrExtWinsRepEnergy.allocate(1);
	ZoneBmSolFrIntWinsRep.allocate(1);
	ZoneBmSolFrIntWinsRepEnergy.allocate(1);

	// Surface
	Surface.allocate(1);
	TotSurfaces = 1;
	SurfaceWindow.allocate(1);
	TotWindows = 1;
	Surface( 1 ).Class = SurfaceClass_Window;
	Surface( 1 ).HeatTransSurf = true;
	// Skip base surface stuff?
	Surface( 1 ).BaseSurf = 1; // It's own base surface?
	Surface( 1 ).ExtBoundCond = 0;
	Surface( 1 ).ExtSolar = true;
	Surface( 1 ).ExtWind = true;
	Surface( 1 ).Zone = 1;
	Surface( 1 ).TAirRef = AdjacentAirTemp;

	SurfaceWindow( 1 ).ShadingFlag = -1;
	SurfaceWindow( 1 ).StormWinFlag = -1;

	HConvIn.allocate(1);
	TempEffBulkAir.allocate(1);
	QHTRadSysSurf.allocate(1);
	QHWBaseboardSurf.allocate(1);
	QSteamBaseboardSurf.allocate(1);
	QElecBaseboardSurf.allocate(1);
	CosIncAng.allocate( 1, 1, 1 );
	SunlitFrac.allocate(1, 1, 1);
	AOSurf.allocate(1);
	SunlitFracWithoutReveal.allocate(1,1,1);
	QRadThermInAbs.allocate(1);
	AirSkyRadSplit.allocate(1);
	QRadSWOutIncident.allocate(1);
	WinHeatGain.allocate(1);
	WinTransSolar.allocate(1);
	WinGainConvGlazToZoneRep.allocate(1);
	WinGainIRGlazToZoneRep.allocate(1);
	WinGapConvHtFlowRep.allocate(1);
	WinGapConvHtFlowRepEnergy.allocate(1);
	QS.allocate(1);
	WinLossSWZoneToOutWinRep.allocate(1);
	WinSysSolTransmittance.allocate(1);
	WinSysSolAbsorptance.allocate(1);
	WinSysSolReflectance.allocate(1);
	InsideGlassCondensationFlag.allocate(1);
	QdotConvOutRep.allocate(1);
	QdotConvOutRepPerArea.allocate(1);
	QConvOutReport.allocate(1);
	QdotRadOutRep.allocate(1);
	QdotRadOutRepPerArea.allocate(1);
	QRadOutReport.allocate(1);
	AISurf.allocate(1);
	AOSurf.allocate(1);
	ISABSF.allocate(1);
	BmIncInsSurfIntensRep.allocate(1);
	BmIncInsSurfAmountRep.allocate(1);
	BmIncInsSurfAmountRepEnergy.allocate(1);
	WinBmSolar.allocate(1);
	WinDifSolar.allocate(1);
	WinBmSolarEnergy.allocate(1);
	WinDifSolarEnergy.allocate(1);
	WinTransSolarEnergy.allocate(1);
	WinBmBmSolar.allocate(1);
	WinBmDifSolar.allocate(1);
	WinBmBmSolarEnergy.allocate(1);
	WinBmDifSolarEnergy.allocate(1);
	WinDirSolTransAtIncAngle.allocate(1);
	AnisoSkyMult.allocate(1);
	CosIncidenceAngle.allocate(1);
	QRadSWOutIncidentBeam.allocate(1);
	QRadSWOutIncidentSkyDiffuse.allocate(1);
	QRadSWOutIncidentGndDiffuse.allocate(1);
	QRadSWOutIncBmToDiffReflGnd.allocate(1);
	QRadSWOutIncSkyDiffReflGnd.allocate(1);
	QRadSWwinAbsTot.allocate(1);
	QRadSWwinAbsTotEnergy.allocate(1);
	SWOutAbsTotalReport.allocate(1);
	SWOutAbsTotalReport.allocate(1);
	WinShadingAbsorbedSolar.allocate(1);

	CosIncAng(1,1,1) = 1.0;
	SunlitFrac(1,1,1) = 1.0;
	SunlitFracWithoutReveal(1,1,1) = 1.0;


}

void remove_dummy_variables()
{
	// Zone
	NumOfZones = 0;
	Zone.deallocate();
	MAT.deallocate();
	ZoneAirHumRatAvg.deallocate();
	ZoneAirHumRat.deallocate();
	DSZone.deallocate();
	DGZone.deallocate();
	DBZoneSSG.deallocate();
	DBZone.deallocate();
	ZoneTransSolar.deallocate();
	ZoneTransSolarEnergy.deallocate();
	ZoneBmSolFrExtWinsRep.deallocate();
	ZoneDifSolFrExtWinsRep.deallocate();
	ZoneBmSolFrExtWinsRepEnergy.deallocate();
	ZoneDifSolFrExtWinsRepEnergy.deallocate();
	ZoneBmSolFrIntWinsRep.deallocate();
	ZoneBmSolFrIntWinsRepEnergy.deallocate();

	// Surface
	Surface.deallocate();
	SurfaceWindow.deallocate();
	TempEffBulkAir.deallocate();
	HConvIn.deallocate();
	QHTRadSysSurf.deallocate();
	QHWBaseboardSurf.deallocate();
	QSteamBaseboardSurf.deallocate();
	QElecBaseboardSurf.deallocate();
	CosIncAng.deallocate();
	SunlitFrac.deallocate();
	AOSurf.deallocate();
	SunlitFracWithoutReveal.deallocate();
	QRadThermInAbs.deallocate();
	AirSkyRadSplit.deallocate();
	QRadSWOutIncident.deallocate();
	WinHeatGain.deallocate();
	WinTransSolar.deallocate();
	WinGainConvGlazToZoneRep.deallocate();
	WinGainIRGlazToZoneRep.deallocate();
	WinGapConvHtFlowRep.deallocate();
	WinGapConvHtFlowRepEnergy.deallocate();
	QS.deallocate();
	WinLossSWZoneToOutWinRep.deallocate();
	WinSysSolTransmittance.deallocate();
	WinSysSolAbsorptance.deallocate();
	WinSysSolReflectance.deallocate();
	InsideGlassCondensationFlag.deallocate();
	QdotConvOutRep.deallocate();
	QdotConvOutRepPerArea.deallocate();
	QConvOutReport.deallocate();
	QdotRadOutRep.deallocate();
	QdotRadOutRepPerArea.deallocate();
	QRadOutReport.deallocate();
	AISurf.deallocate();
	AOSurf.deallocate();
	ISABSF.deallocate();
	BmIncInsSurfIntensRep.deallocate();
	BmIncInsSurfAmountRep.deallocate();
	BmIncInsSurfAmountRepEnergy.deallocate();
	WinBmSolar.deallocate();
	WinDifSolar.deallocate();
	WinBmSolarEnergy.deallocate();
	WinDifSolarEnergy.deallocate();
	WinTransSolarEnergy.deallocate();
	WinBmBmSolar.deallocate();
	WinBmDifSolar.deallocate();
	WinBmBmSolarEnergy.deallocate();
	WinBmDifSolarEnergy.deallocate();
	WinDirSolTransAtIncAngle.deallocate();
	AnisoSkyMult.deallocate();
	CosIncidenceAngle.deallocate();
	QRadSWOutIncidentBeam.deallocate();
	QRadSWOutIncidentSkyDiffuse.deallocate();
	QRadSWOutIncidentGndDiffuse.deallocate();
	QRadSWOutIncBmToDiffReflGnd.deallocate();
	QRadSWOutIncSkyDiffReflGnd.deallocate();
	QRadSWwinAbsTot.deallocate();
	QRadSWwinAbsTotEnergy.deallocate();
	SWOutAbsTotalReport.deallocate();
	SWOutAbsTotalReport.deallocate();
	WinShadingAbsorbedSolar.deallocate();

	// Environment
	BeamSolarRad = 0.0;
	SunIsUp = false;

}

} // WindowASHRAE1588RP

} // EnergyPlus
