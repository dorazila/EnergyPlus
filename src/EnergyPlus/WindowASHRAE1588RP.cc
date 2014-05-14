// C++ Headers
#include <string>

// ObjexxFCL Headers
#include <ObjexxFCL/gio.hh>

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
#include <DataSystemVariables.hh>
#include <DisplayRoutines.hh>
#include <DataTimings.hh>
#include <General.hh>
#include <HeatBalanceManager.hh>
#include <HeatBalanceSurfaceManager.hh>
#include <InputProcessor.hh>
#include <WindowManager.hh>
#include <SolarShading.hh>
#include <UtilityRoutines.hh>

namespace EnergyPlus {

namespace WindowASHRAE1588RP {

using namespace DataEnvironment;
using namespace DataGlobals;
using namespace DataHeatBalance;
using namespace DataHeatBalFanSys;
using namespace DataHeatBalSurface;
using namespace DataIPShortCuts;
using namespace DataSurfaces;
using namespace DataSystemVariables;
using namespace DataTimings;
using namespace General;

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
	FArray1D_string ConstructAlphas( 7 ); // Construction Alpha names defined
	FArray1D< Real64 > ConstructNumerics( 6 ); // Temporary array to transfer construction properties
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

		int TotMaterialsSave = TotMaterials;

		// Save Materials
		MaterialSave.allocate( TotMaterials );
		NominalRSave.allocate( TotMaterials );
		MaterialSave = Material;
		NominalRSave = NominalR;
		Material.deallocate();
		NominalR.deallocate();

		int number_of_gaps;
		int number_of_new_materials;

		FArray1D< MaterialProperties > new_materials;
		FArray1D< Real64 > new_nominal_R;


		FArray1D< Real64 > NominalRforNominalUCalculationSave;
		int TotConstructsSave = TotConstructs;

		// Save Constructions -- The list will be deleted so that the only
		// construction is the one currently being set for any borrowed subroutines
		{
			ConstructSave.allocate( TotConstructs );
			ConstructSave = Construct;
			NominalRforNominalUCalculationSave.allocate( TotConstructs );
			NominalUSave.allocate( TotConstructs );
			NominalRforNominalUCalculationSave = NominalRforNominalUCalculation;
			NominalUSave = NominalU;

			Construct.deallocate();
			NominalRforNominalUCalculation.deallocate();
			NominalU.deallocate();

			Construct.allocate(1);
			NominalRforNominalUCalculation.allocate(1);
			NominalU.allocate(1);

			TotConstructs = 1;
		}


		ConstructionData new_construct;

		new_construct.Name = ConstructAlphas( 1 );
		new_construct.TypeIsWindow = true;

		// Save Frame and Divider objects
		int TotFrameDividerSave = TotFrameDivider;
		FArray1D< FrameDividerProperties > FrameDividerSave;

		FrameDividerSave.allocate( TotFrameDivider );
		FrameDividerSave = FrameDivider;

		FrameDivider.deallocate();

		FrameDivider.allocate(1);

		TotFrameDivider = 1;

		FrameDividerProperties new_frame_divider;


		// detect analysis type (stand-alone vs. integrated)
		bool stand_alone_analysis;
		std::string ashrae1588_file_name;

		if ( lAlphaFieldBlanks( 7 ) )
		{
			stand_alone_analysis = false;
			ashrae1588_file_name = "";
		}
		else
		{
			stand_alone_analysis = true;
			ashrae1588_file_name = ConstructAlphas(7);
		}


		// set initial values and set locks as appropriate. If IDF is blank,
		// override with defaults from ASHRAE 1588 RP Database file (TODO)
		Real64 target_u_factor;
		bool u_factor_set;

		if ( lNumericFieldBlanks( 1 ) )
		{
			u_factor_set = false;
		}
		else
		{
			u_factor_set = true;
			target_u_factor = ConstructNumerics( 1 );
		}

		Real64 target_shgc;
		bool shgc_set;

		if ( lNumericFieldBlanks( 2 ) )
		{
			shgc_set = false;
		}
		else
		{
			shgc_set = true;
			target_shgc = ConstructNumerics( 2 );
		}

		Real64 target_vt;
		bool vt_set;

		if ( lNumericFieldBlanks( 3 ) )
		{
			vt_set = false;
		}
		else
		{
			vt_set = true;
			target_vt = ConstructNumerics( 3 );
		}

		std::string fenestration_type;
		bool fenestration_type_lock;

		if ( lAlphaFieldBlanks( 2 ) )
		{
			fenestration_type_lock = false;
			fenestration_type = "HORIZONTALSLIDER";
		}
		else
		{
			fenestration_type_lock = true;
			fenestration_type = ConstructAlphas( 2 );
		}

		int number_of_panes;
		bool number_of_panes_lock;

		if ( lNumericFieldBlanks( 4 ) )
		{
			number_of_panes_lock = false;
			number_of_panes = 2;
		}
		else
		{
			number_of_panes_lock = true;
			number_of_panes = ConstructNumerics( 4 );
		}

		std::string glazing_type;
		bool glazing_type_lock;

		if ( lAlphaFieldBlanks( 3 ) )
		{
			glazing_type_lock = false;
			glazing_type = "CLEAR";
		}
		else
		{
			glazing_type_lock = true;
			glazing_type = ConstructAlphas( 3 );
		}

		std::string glazing_surface_treatment;
		bool glazing_surface_treatment_lock;

		if ( lAlphaFieldBlanks( 4 ) )
		{
			glazing_surface_treatment_lock = false;
			glazing_surface_treatment = "NONE";
		}
		else
		{
			glazing_surface_treatment_lock = true;
			glazing_surface_treatment = ConstructAlphas( 4 );
		}

		std::string gas_type;
		bool gas_type_lock;

		if ( lAlphaFieldBlanks( 5 ) )
		{
			gas_type_lock = false;
			gas_type = "AIR";
		}
		else
		{
			gas_type_lock = true;
			gas_type = ConstructAlphas( 5 );
		}

		std::string frame_material;
		bool frame_material_lock;

		if ( lAlphaFieldBlanks( 6 ) )
		{
			frame_material_lock = false;
			frame_material = "VINYL";
		}
		else
		{
			frame_material_lock = true;
			frame_material = ConstructAlphas( 6 );
		}

		Real64 frame_width;
		bool frame_width_lock;

		if ( lNumericFieldBlanks( 5 ) )
		{
			frame_width_lock = false;
			frame_width = 0.0;
		}
		else
		{
			frame_width_lock = true;
			frame_width = ConstructNumerics( 5 );
		}

		Real64 divider_width;
		bool divider_width_lock;

		if ( lNumericFieldBlanks( 6 ) )
		{
			divider_width_lock = false;
			divider_width = 0.0;
		}
		else
		{
			divider_width_lock = true;
			divider_width = ConstructNumerics( 6 );
		}



		// internal defaults to be varied. TODO read these from ASHRAE 1588 RP Database file, or derive them as appropriate from other inputs.
		Real64 glass_thickness = 0.003;
		Real64 glass_solar_transmissivity = 0.837;
		Real64 glass_visible_transmissivity = 0.898;
		Real64 glass_solar_reflectivity = 0.075;
		Real64 glass_visible_reflectivity = 0.081;
		Real64 glass_IR_transmissivity = 0.0;
		Real64 glass_IR_absorptivity = 0.84;

		Real64 gap_thickness = 0.0127;

		Real64 frame_solar_absorptivity = 0.7;
		Real64 frame_visible_absorptivity = 0.7;
		Real64 frame_IR_emissivity = 0.7;


		// internal defaults based on other values
		Real64 frame_conductance = 10.0;
		Real64 frame_edge_ratio = 1.0;

		Real64 fenestration_width;
		Real64 fenestration_height;
		Real64 glazing_width;
		Real64 glazing_height;
		Real64 tilt;
		Real64 fenestration_area;
		Real64 glazing_area;

		Real64 glass_conductivity;

		bool has_frame;
		int num_horizontal_dividers;
		int num_vertical_dividers;

		// matching variables
		Real64 u_factor;
		Real64 shgc;
		Real64 vt;

		// internal defaults to be left alone
		Real64 glass_youngs_modulus = 7.2e10;
		Real64 glass_poissons_ratio = 0.22;

		Real64 max_divider_spacing = 0.3; // NFRC 100-2014 4.2.2 (B)

		// Allocate temporary arrays
		create_dummy_variables();

		Surface( 1 ).Name = ConstructAlphas( 1 ) + ":Surface";

		ASHRAE1588RP_Flag = true;
		KickOffSimulation = false;

		bool target_matched = false;
		// This is where the iterative optimization loop will begin
		while (! target_matched)
		{

			// internal defaults based on other values

			// set product sizes and tilts based on NFRC 100-2014 Table 4-3
			if ( fenestration_type == "CASEMENTDOUBLE" )
			{
				fenestration_width = 1.2;
				fenestration_height = 1.5;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "CASEMENTSINGLE" )
			{
				fenestration_width = 0.6;
				fenestration_height = 1.5;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "DUALACTION" )
			{
				fenestration_width = 1.2;
				fenestration_height = 1.5;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "FIXED" )
			{
				fenestration_width = 1.2;
				fenestration_height = 1.5;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "GARAGEORROLLINGDOOR" )
			{
				fenestration_width = 2.134;
				fenestration_height = 2.134;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "GREENHOUSEORGARDEN" )
			{
				fenestration_width = 1.5;
				fenestration_height = 1.2;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "HINGEDESCAPE" )
			{
				fenestration_width = 1.5;
				fenestration_height = 1.2;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "HORIZONTALSLIDER" )
			{
				fenestration_width = 1.5;
				fenestration_height = 1.2;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "HYBRIDTUBULARDAYLIGHTINGDEVICE" )
			{
				fenestration_width = 0.4697;
				fenestration_height = 0.4697;
				tilt = 0.0; // 0 deg
			}
			else if ( fenestration_type == "JALORJALAWNING" )
			{
				fenestration_width = 1.2;
				fenestration_height = 1.5;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "PIVOTED" )
			{
				fenestration_width = 1.2;
				fenestration_height = 1.5;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "PROJECTINGAWNINGDUAL" )
			{
				fenestration_width = 1.5;
				fenestration_height = 1.2;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "PROJECTINGAWNINGSINGLE" )
			{
				fenestration_width = 1.5;
				fenestration_height = 0.6;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "DOORSIDELITE" )
			{
				fenestration_width = 0.6;
				fenestration_height = 2.0;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "SKYLIGHTORROOFWINDOW" )
			{
				fenestration_width = 1.2;
				fenestration_height = 1.2;
				tilt = Pi/9; // 20 deg
			}
			else if ( fenestration_type == "SLIDINGPATIODOORWITHFRAME" )
			{
				fenestration_width = 2.0;
				fenestration_height = 2.0;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "CURTAINWALLORWINDOWWALL" )
			{
				fenestration_width = 2.0;
				fenestration_height = 2.0;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "SLOPEDGLAZING" )
			{
				fenestration_width = 2.0;
				fenestration_height = 2.0;
				tilt = Pi/9; // 20 deg
			}
			else if ( fenestration_type == "SPANDRELPANEL" )
			{
				fenestration_width = 2.0;
				fenestration_height = 1.2;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "SWINGINGDOORWITHFRAME" )
			{
				// Assume single door (double door width = 1.92 m)
				fenestration_width = 0.96;
				fenestration_height = 2.09;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "DOORTRANSOM" )
			{
				fenestration_width = 2.0;
				fenestration_height = 0.6;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "TROPICALAWNING" )
			{
				fenestration_width = 1.5;
				fenestration_height = 1.2;
				tilt = Pi/2; // 90 deg
			}
			else if ( fenestration_type == "TUBULARDAYLIGHTINGDEVICE" )
			{
				fenestration_width = 0.3102;
				fenestration_height = 0.3102;
				tilt = 0.0; // 0 deg
			}
			else if ( fenestration_type == "VERTICALSLIDER" )
			{
				fenestration_width = 1.2;
				fenestration_height = 1.5;
				tilt = Pi/2; // 90 deg
			}

			if ( glazing_type == "POLYESTERFILM" )
			{
				glass_conductivity = 0.14;
			}
			else
			{
				glass_conductivity = 0.9;
			}

			if ( frame_width > 0.0 )
			{
				has_frame = true;
			}
			else
			{
				has_frame = false;
			}

			fenestration_area = fenestration_width*fenestration_height;
			glazing_width = fenestration_width - 2.0*frame_width;
			glazing_height = fenestration_height - 2.0*frame_width;
			glazing_area = glazing_width*glazing_height;

			if ( has_frame )
			{
				num_horizontal_dividers = ceil(glazing_height/max_divider_spacing);
				num_vertical_dividers = ceil(glazing_width/max_divider_spacing);
			}
			else
			{
				num_horizontal_dividers = 0;
				num_vertical_dividers = 0;
			}


			Surface( 1 ).Height = glazing_height;
			Surface( 1 ).Width = glazing_width;
			Surface( 1 ).Area = glazing_area;
			Surface( 1 ).Tilt = tilt*180/Pi;
			Surface( 1 ).CosTilt = cos(tilt);
			Surface( 1 ).SinTilt = sin(tilt);
			Surface( 1 ).ViewFactorSky = 0.5 * ( 1.0 + Surface( 1 ).CosTilt );
			Surface( 1 ).ViewFactorGround = 0.5 * ( 1.0 - Surface( 1 ).CosTilt );
			Surface( 1 ).ViewFactorSkyIR = Surface( 1 ).ViewFactorSky;
			Surface( 1 ).ViewFactorGroundIR = Surface( 1 ).ViewFactorGround;
			AirSkyRadSplit( 1 ) = std::sqrt( 0.5 * ( 1.0 + Surface( 1 ).CosTilt ) );

			number_of_gaps = number_of_panes - 1;
			number_of_new_materials = number_of_panes + number_of_gaps;

			// Construction specific allocations
			AWinSurf.allocate(1, number_of_panes);
			QRadSWwinAbs.allocate(1, number_of_panes);
			QRadSWwinAbsLayer.allocate(1, number_of_panes);

			// Create New Material objects
			if ( new_materials.size_ != (unsigned)number_of_new_materials ) {
				new_materials.allocate( number_of_new_materials );
				Material.allocate( number_of_new_materials );
				NominalR.allocate( number_of_new_materials );
				TotMaterials = number_of_new_materials;
			}


			// Define material properties for glazings
			for ( int MaterNum = 1; MaterNum <= number_of_new_materials; MaterNum += 2 )
			{
				Material( MaterNum ).Group = WindowGlass;
				Material( MaterNum ).Name = ConstructAlphas( 1 ) + ":GLAZING" + std::to_string(MaterNum);
				Material( MaterNum ).Roughness = VerySmooth;
				Material( MaterNum ).ROnly = true;
				Material( MaterNum ).Thickness = glass_thickness;
				Material( MaterNum ).Trans = glass_solar_transmissivity;
			    Material( MaterNum ).ReflectSolBeamFront = glass_solar_reflectivity;
				Material( MaterNum ).ReflectSolBeamBack = glass_solar_reflectivity;
				Material( MaterNum ).TransVis = glass_visible_transmissivity;
				Material( MaterNum ).ReflectVisBeamFront = glass_visible_reflectivity;
				Material( MaterNum ).ReflectVisBeamBack = glass_visible_reflectivity;
				Material( MaterNum ).TransThermal = glass_IR_transmissivity;
				Material( MaterNum ).AbsorpThermalFront = glass_IR_absorptivity;
				Material( MaterNum ).AbsorpThermalBack = glass_IR_absorptivity;
				Material( MaterNum ).Conductivity = glass_conductivity;
				Material( MaterNum ).GlassTransDirtFactor = 1.0;  // TODO Expose?
				Material( MaterNum ).YoungModulus = glass_youngs_modulus;
				Material( MaterNum ).PoissonsRatio = glass_poissons_ratio;
				Material( MaterNum ).AbsorpThermal = Material( MaterNum ).AbsorpThermalBack;
				Material( MaterNum ).SolarDiffusing = false;  // TODO Expose?

				Material( MaterNum ).GlassSpectralDataPtr = 0;

				NominalR( MaterNum ) = Material( MaterNum ).Thickness / Material( MaterNum ).Conductivity;
				Material( MaterNum ).Resistance = NominalR( MaterNum );

			}

			// Define material properties for gaps
			for ( int MaterNum = 2; MaterNum <= number_of_new_materials; MaterNum += 2 )
			{
				Material( MaterNum ).Group = WindowGas;
				Material( MaterNum ).Name = ConstructAlphas( 1 ) + ":GAP" + std::to_string(MaterNum);
				Material( MaterNum ).Roughness = MediumRough;
				Material( MaterNum ).ROnly = true;
				Material( MaterNum ).Thickness = gap_thickness;
				Material( MaterNum ).NumberOfGasesInMixture = 1;
				Material( MaterNum ).GasFract( 1 ) = 1.0;


				if ( gas_type == "AIR" ) Material( MaterNum ).GasType( 1 ) = 1;
				if ( gas_type == "ARGON" ) Material( MaterNum ).GasType( 1 ) = 2;
				if ( gas_type == "KRYPTON" ) Material( MaterNum ).GasType( 1 ) = 3;
				if ( gas_type == "XENON" ) Material( MaterNum ).GasType( 1 ) = 4;

				Material( MaterNum ).GasWght( 1 ) = GasWght( Material( MaterNum ).GasType( 1 ) );
				Material( MaterNum ).GasSpecHeatRatio( 1 ) = GasSpecificHeatRatio( Material( MaterNum ).GasType( 1 ) );
				for ( int ICoeff = 1; ICoeff <= 3; ++ICoeff ) {
					Material( MaterNum ).GasCon( 1, ICoeff ) = GasCoeffsCon( Material( MaterNum ).GasType( 1 ), ICoeff );
					Material( MaterNum ).GasVis( 1, ICoeff ) = GasCoeffsVis( Material( MaterNum ).GasType( 1 ), ICoeff );
					Material( MaterNum ).GasCp( 1, ICoeff ) = GasCoeffsCp( Material( MaterNum ).GasType( 1 ), ICoeff );
				}

				Real64 DenomRGas = ( Material( MaterNum ).GasCon( 1, 1 ) + Material( MaterNum ).GasCon( 1, 2 ) * 300.0 + Material( MaterNum ).GasCon( 1, 3 ) * 90000.0 );
				NominalR( MaterNum ) = Material( MaterNum ).Thickness / DenomRGas;

			}

			new_construct.TotLayers = number_of_new_materials;

			for ( int Layer = 1; Layer <= number_of_new_materials; ++Layer ) {
				new_construct.LayerPoint( Layer ) = Layer;
			}

			Construct( 1 ) = new_construct;

			NominalRforNominalUCalculation( 1 ) = 0.0;
			for ( int Layer = 1; Layer <= Construct( 1 ).TotLayers; ++Layer ) {
				NominalRforNominalUCalculation( 1 ) += NominalR( Construct( 1 ).LayerPoint( Layer ) );
			}

			CheckAndSetConstructionProperties( 1, ErrorsFound );

			// Set frame and divider properties
			if ( has_frame )
			{
				new_frame_divider.Name = ConstructAlphas( 1 ) + ":FRAME";
				new_frame_divider.FrameWidth = frame_width;
				new_frame_divider.FrameProjectionOut = 0.0;
				new_frame_divider.FrameProjectionIn = 0.0;
				new_frame_divider.FrameConductance = frame_conductance;
				new_frame_divider.FrEdgeToCenterGlCondRatio = frame_edge_ratio;
				new_frame_divider.FrameSolAbsorp = frame_solar_absorptivity;
				new_frame_divider.FrameVisAbsorp = frame_visible_absorptivity;
				new_frame_divider.FrameEmis = frame_IR_emissivity;
				new_frame_divider.FrameEdgeWidth = 0.06355; // 2.5 in
				new_frame_divider.DividerType = DividedLite;
				new_frame_divider.DividerWidth = divider_width;
				new_frame_divider.HorDividers = num_horizontal_dividers;
				new_frame_divider.VertDividers = num_vertical_dividers;
				new_frame_divider.DividerProjectionOut = 0.0;
				new_frame_divider.DividerProjectionIn = 0.0;
				new_frame_divider.DividerConductance = frame_conductance;
				new_frame_divider.DivEdgeToCenterGlCondRatio = frame_edge_ratio;
				new_frame_divider.DividerSolAbsorp = frame_solar_absorptivity;
				new_frame_divider.DividerVisAbsorp = frame_visible_absorptivity;
				new_frame_divider.DividerEmis = frame_IR_emissivity;
				new_frame_divider.DividerEdgeWidth = 0.06355; // 2.5 in

				SurfaceWindow( 1 ).FrameArea = fenestration_area - glazing_area;
				SurfaceWindow( 1 ).DividerArea = divider_width*(num_horizontal_dividers*glazing_width + num_vertical_dividers*glazing_height - num_horizontal_dividers*num_vertical_dividers*divider_width);
				Surface( 1 ).Area -= SurfaceWindow( 1 ).DividerArea;
				SurfaceWindow( 1 ).GlazedFrac = Surface( 1 ).Area / ( Surface( 1 ).Area + SurfaceWindow( 1 ).DividerArea );


				FrameDivider( 1 ) = new_frame_divider;

				Surface( 1 ).FrameDivider = 1;
			}
			else
			{
				Surface( 1 ).FrameDivider = 0;
			}

			Surface( 1 ).Construction = 1; // This is the only construction available to the dummy surface. The actual surface will reference the real construction.


			// Setup functions

			InitGlassOpticalCalculations();

			// Set up U-factor conditions
			Real64 in_air_temp = 21.0;
			Real64 out_air_temp = -18.0;
			Real64 wind_speed = 5.5;
			Real64 solar_inccident = 0.0;

			calc_window_performance(in_air_temp, out_air_temp, wind_speed, solar_inccident);

			u_factor = -WinHeatGain(1)/(fenestration_area*(in_air_temp - out_air_temp));

			// Set up SHGC conditions
			in_air_temp = 24.0;
			out_air_temp = 32.0;
			wind_speed = 2.75;
			solar_inccident = 783.0;

			calc_window_performance(in_air_temp, out_air_temp, wind_speed, solar_inccident);

			Real64 q_total = WinHeatGain(1);

			// NFRC 201-2014 Equation 8-7
			Real64 q_U = u_factor*fenestration_area*(out_air_temp - in_air_temp);

			// NFRC 201-2014 Equation 8-2
			shgc = (q_total - q_U)/(fenestration_area*solar_inccident);

			Real64 non_opaque_area_fraction = Surface( 1 ).Area/fenestration_area;
			vt = POLYF(1.0,Construct( 1 ).TransVisBeamCoef( 1 ))*non_opaque_area_fraction;

			// if match not obtained adjust inputs

			// Deallocate construction specific arrays
			AWinSurf.deallocate();
			QRadSWwinAbs.deallocate();
			QRadSWwinAbsLayer.deallocate();

			if (!u_factor_set && !shgc_set && !vt_set) target_matched = true;

			target_matched = true;

		} // matching loop

		ASHRAE1588RP_Flag = false;
		KickOffSimulation = true;

		// TODO This will trigger on the first ASHRAE 1588 construction found.
		// It should be moved out of the construction loop to allow for multiple
		// constructions.
		if ( stand_alone_analysis )
		{
			// Write to file
			int output_file;
			output_file = GetNewUnitNumber();

			{ IOFlags flags; flags.ACTION( "write" ); gio::open( output_file, ashrae1588_file_name, flags );}
			gio::write( output_file, "(A)" ) << "Target U-factor: " + RoundSigDigits(target_u_factor,3);
			gio::write( output_file, "(A)" ) << "Target SHGC: " + RoundSigDigits(target_shgc,3);
			gio::write( output_file, "(A)" ) << "Target Visible Transmittance: " + RoundSigDigits(target_vt,3);

			gio::write( output_file, "(A)" ) << "\nMatch U-factor: " + RoundSigDigits(u_factor,3);
			gio::write( output_file, "(A)" ) << "Match SHGC: " + RoundSigDigits(shgc,3);
			gio::write( output_file, "(A)" ) << "Match Visible Transmittance: " + RoundSigDigits(vt,3);

			gio::close( output_file );



			// Write to console
			std::string Elapsed;
			int Hours; // Elapsed Time Hour Reporting
			int Minutes; // Elapsed Time Minute Reporting
			Real64 Seconds; // Elapsed Time Second Reporting
			gio::Fmt ETimeFmt( "(I2.2,'hr ',I2.2,'min ',F5.2,'sec')" );

			Time_Finish = epElapsedTime();
			if ( Time_Finish < Time_Start ) Time_Finish += 24.0 * 3600.0;
			Elapsed_Time = Time_Finish - Time_Start;
			Hours = Elapsed_Time / 3600.0;
			Elapsed_Time -= Hours * 3600.0;
			Minutes = Elapsed_Time / 60.0;
			Elapsed_Time -= Minutes * 60.0;
			Seconds = Elapsed_Time;
			if ( Seconds < 0.0 ) Seconds = 0.0;
			gio::write( Elapsed, ETimeFmt ) << Hours << Minutes << Seconds;

			gio::write( "(1X,A)" ) << ( "EnergyPlus ASHRAE 1588-RP Window Construction Generated Successfully-- Elapsed Time=" + Elapsed );
			exit (EXIT_SUCCESS);
		}

		// deallocate temporary arrays
		remove_dummy_variables();

		// Restore materials list and copy in new materials
		{
			new_materials = Material;
			new_nominal_R = NominalR;

			Material.deallocate();
			NominalR.deallocate();

			TotMaterials = TotMaterialsSave;

			Material.allocate( TotMaterials + number_of_new_materials);
			NominalR.allocate( TotMaterials + number_of_new_materials);
			Material( {1,TotMaterials} ) = MaterialSave( {1,TotMaterials} );
			NominalR( {1,TotMaterials} ) = NominalRSave( {1,TotMaterials} );
			Material( {TotMaterials + 1, TotMaterials + number_of_new_materials} ) = new_materials;
			NominalR( {TotMaterials + 1, TotMaterials + number_of_new_materials} ) = new_nominal_R;

			MaterialSave.deallocate();
			NominalRSave.deallocate();
		}

		// Restore frame and divider list and copy in new frame and divider
		FrameDivider.deallocate();

		TotFrameDivider = TotFrameDividerSave;

		if ( has_frame )
		{
			FrameDivider.allocate( TotFrameDivider + 1 );
			FrameDivider( {1,TotFrameDivider} ) = FrameDividerSave;
			FrameDivider( TotFrameDivider + 1 ) = new_frame_divider;
			TotFrameDivider += 1;
		}
		else
		{
			FrameDivider.allocate( TotFrameDivider );
			FrameDivider = FrameDividerSave;
		}

		FrameDividerSave.deallocate();


		// Restore construction list and copy in new construction
		{
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
			NominalU = NominalUSave;


			Construct( ConstrNum ) = new_construct;
			// Set new layer references corresponding to new material numbers
			for ( int Layer = 1; Layer <= Construct( ConstrNum ).TotLayers; ++Layer ) {
				Construct( ConstrNum ).LayerPoint( Layer ) = TotMaterials + Layer;
			}

			TotMaterials += number_of_new_materials;

			if ( has_frame )
			{
				Construct( ConstrNum ).W5FrameDivider = TotFrameDivider;
			}

			NominalRforNominalUCalculation( ConstrNum ) = newR;
			NominalU( ConstrNum ) = newU;

			ConstructSave.deallocate();
			NominalRforNominalUCalculationSave.deallocate();
			NominalUSave.deallocate();
		}


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

	// Save tilt information for natural convection calculations
	Real64 tilt_save = Surface( 1 ).Tilt;

	for (int i = 0; i < max_iterations; i++) {

		// Use complementary angle for exterior natural convection calculations
		Surface( 1 ).Tilt = 180 - tilt_save;
		Surface( 1 ).CosTilt = cos(Surface( 1 ).Tilt*Pi/180);
		Surface( 1 ).SinTilt = sin(Surface( 1 ).Tilt*Pi/180);
		CalcISO15099WindowIntConvCoeff( 1, out_surf_temp, T_out); // This subroutine sets the global HConvIn( 1 ) variable. We will use it to set the exterior natural convection.
		h_exterior = h_exterior_f + HConvIn( 1 ); // add natural convection

		// revert tilt for interior natural convection calculations
		Surface( 1 ).Tilt = tilt_save;
		Surface( 1 ).CosTilt = cos(tilt_save*Pi/180);
		Surface( 1 ).SinTilt = sin(tilt_save*Pi/180);
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

	ZoneAirHumRatAvg(1) = 0.0;
	ZoneAirHumRat(1) = 0.0;

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
	WinGainFrameDividerToZoneRep.allocate(1);
	InsideFrameCondensationFlag.allocate(1);
	InsideDividerCondensationFlag.allocate(1);

	QHTRadSysSurf(1) = 0.0;
	QHWBaseboardSurf(1) = 0.0;
	QSteamBaseboardSurf(1) = 0.0;
	QElecBaseboardSurf(1) = 0.0;
	QRadThermInAbs(1) = 0.0;
	QS(1) = 0.0;
	ISABSF(1) = 0.0;
	CosIncAng(1,1,1) = 1.0;
	SunlitFrac(1,1,1) = 1.0;
	SunlitFracWithoutReveal(1,1,1) = 1.0;
	AnisoSkyMult(1) = 0.0;  // This may need to change if NFRC adds a diffuse component for SHGC tests



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
	WinGainFrameDividerToZoneRep.deallocate();
	InsideFrameCondensationFlag.deallocate();
	InsideDividerCondensationFlag.deallocate();

	// Environment
	BeamSolarRad = 0.0;
	SunIsUp = false;

}

} // WindowASHRAE1588RP

} // EnergyPlus
