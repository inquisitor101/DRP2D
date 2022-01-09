#include "config_structure.hpp"



CConfig::CConfig
(
 const char *configFile
)
 /*
	* Constructor, reads the specified input configuration file.
	*/
{
  // Message stream.
	std::ostringstream message;

  // Check if file exists.
  std::ifstream inputFile(configFile);
  if( !inputFile.good() )
    Terminate("CConfig::CConfig", __FILE__, __LINE__,
              "File could not be opened!");


	if( !ReadGridOptions(configFile) ){
		message << "Failed to extract grid options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

 	// Extract input/output information.
	if( !ReadIOOptions(configFile) ){
		message << "Failed to read input/output options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

	// Extract boundary marker specification.
	if( !ReadBoundaryOptions(configFile) ){
		message << "Failed to read boundary marker options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract temporal information.
	if( !ReadTemporalOptions(configFile) ){
		message << "Failed to read temporal options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract flow information.
	if( !ReadFlowOptions(configFile) ){
		message << "Failed to read flow options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract initial condition information.
	if( !ReadICOptions(configFile) ){
		message << "Failed to read initial condition options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

  // Extract solver options.
  if( !ReadSolverOptions(configFile) ){
		message << "Failed to read solver options from " << configFile;
		Terminate("CConfig::CConfig", __FILE__, __LINE__, message.str());
	}

	// Close file.
  inputFile.close();
}


CConfig::~CConfig
(
 void
)
 /*
	* Destructor for CConfig, does nothing.
	*/
{

}


bool CConfig::ReadSolverOptions
(
 const char *configFile
)
 /*
	* Function that reads the solver specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read type of stencils used in first-derivative approximation.
	AddVectorOption(paramFile, "TYPE_STENCIL", NameTypeStencil, true);
  // Assign TypeStencil map.
  MapTypeStencil();

  // Read buffer layer type.
  AddVectorOption(paramFile, "TYPE_BUFFER_LAYER", NameTypeBufferLayer, DefaultParam.NameTypeBufferLayer, true);
  // Assign BufferLayerType map.
  MapTypeBufferLayer();

  // Consistency check.
  if( NameTypeStencil.size() != nDim )
    Terminate("CConfig::ReadSolverOptions", __FILE__, __LINE__,
              "TYPE_STENCIL must be of dimension: 2.");


  // Close file.
  paramFile.close();

	// Return happily.
	return true;
}


bool CConfig::ReadGridOptions
(
 const char *configFile
)
 /*
	* Function that reads all grid information.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Number of zones expected.
	AddScalarOption(paramFile, "NUMBER_ZONE", nZone, DefaultParam.nZone, true);

  // Read zone regions.
	AddVectorOption(paramFile, "MARKER_ZONE", NameZoneMarker, DefaultParam.NameZoneMarker, true);
  // Pad entries for NameZoneMarker.
  PadEntriesVectorData(NameZoneMarker, "MARKER_ZONE", nZone);

  // Read domain bounds in each zone.
	AddVectorOption(paramFile, "DOMAIN_BOUND", DomainBound, true);

	// Read number of nodes in x-direction, per zone.
	AddVectorOption(paramFile, "NUMBER_XNODE", nxNode, true);
  // Pad entries for nxNode.
  PadEntriesVectorData(nxNode, "NUMBER_XNODE", nZone, 1, 2);

	// Read number of nodes in y-direction, per zone.
	AddVectorOption(paramFile, "NUMBER_YNODE", nyNode, true);
  // Pad entries for nyNode.
  PadEntriesVectorData(nyNode, "NUMBER_YNODE", nZone, 1, 2);

  // Read zone conformity, if specified.
  AddBoolOption(paramFile, "ZONE_CONFORMITY", ZoneConformity, DefaultParam.NameZoneConformity, true);
  // Read nodal ratio sizes.
  AddVectorOption(paramFile, "NODAL_RATIO", dsNodalRatioZone, DefaultParam.dsNodalRatioZone, true);
  // Read whether a uniform grid is used or not.
  AddBoolOption(paramFile, "UNIFORM_GRID_RESOLUTION", UniformGridResolution, DefaultParam.NameUniformGridResolution, true);


  // Assign ZoneMarker.
  MapTypeZone();
  // Determine multizone strategy to use.
  DetermineMultizoneStrategy();
  // Check nodal ratio specified according to multizone strategy.
  CheckdsNodalRatio();

  // Process zone conformity option so the solver overwrites input number of
  // nodes in the different zones so the grid remains consistent with ZONE_MAIN.
  if( ZoneConformity ) ProcessZoneConformity();

  // Assign solver type. For now, use a EE only.
  TypeSolver.resize(nZone);
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    TypeSolver[iZone] = SOLVER_EE;

   // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


bool CConfig::ReadICOptions
(
 const char *configFile
)
 /*
	* Function that reads the initial conditions.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

	// Read initial condition in regions.
	AddScalarOption(paramFile, "TYPE_IC", NameInitialCondition, true);
	// Assign TypeIC map.
	MapTypeIC();

  // Read the disturbance center, if specified.
  AddVectorOption(paramFile, "DISTURBANCE_CENTER", CenterX0, DefaultParam.CenterX0, true);

  // Read the disturbance peak percentage, if specified.
  AddScalarOption(paramFile, "DISTURBANCE_RATIO", DisturbanceRatio, DefaultParam.DisturbanceRatio, true);
  // Read the disturbance width, if specified.
  AddScalarOption(paramFile, "DISTURBANCE_WIDTH", DisturbanceWidth, DefaultParam.DisturbanceWidth, true);

	// Check disturbance center, must be a multiple of nDim.
	if( CenterX0.size()%2 != 0)
	  Terminate("CConfig::ReadICOptions", __FILE__, __LINE__,
	            "CenterX0 must be of size multiple of nDim,");

  // Close file.
  paramFile.close();

	// Return happily.
	return true;
}


bool CConfig::ReadFlowOptions
(
 const char *configFile
)
 /*
	* Function that reads the flow conditions.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read the free-stream Mach number, if specified.
  AddScalarOption(paramFile, "FREESTREAM_MACH", MachInf, DefaultParam.MachInf, true);
  // Read the flow angle, if specified.
  AddScalarOption(paramFile, "FLOW_ANGLE", FlowAngle, DefaultParam.FlowAngle, true);

  // Check whether a cross-flow is present of not.
  CrossFlow = ( fabs(sin(FlowAngle*PI_CONSTANT/180.0)) < 1e-10 ) ? false : true;

  // Check if cross-flow is present.
  if( CrossFlow )
    Terminate("CConfig::ReadFlowOptions", __FILE__, __LINE__,
              "Cross-flow is not currently supported (yet).");

  // Close file.
  paramFile.close();

	// Return happily.
	return true;
}


bool CConfig::ReadIOOptions
(
 const char *configFile
)
 /*
	* Function that reads the output information.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read restart solution, if specified.
  AddBoolOption(paramFile, "RESTART_SOLUTION", RestartSolution, DefaultParam.NameRestartSolution, true);
  // Read output file-writing frequency.
  AddScalarOption(paramFile, "WRITE_FREQ", WriteFreq, DefaultParam.WriteFreq, true);
  // Read output screen-monitoring frequency.
  AddScalarOption(paramFile, "OUTPUT_FREQ", OutputFreq, DefaultParam.OutputFreq, true);
	// Read output solution filename.
	AddScalarOption(paramFile, "OUTPUT_SOL_FILENAME", OutputSolFilename, true);
	// Read output solution filename.
	AddScalarOption(paramFile, "OUTPUT_VTK_FILENAME", OutputVTKFilename, true);
  // Read input solution restart filename.
	AddScalarOption(paramFile, "RESTART_FILENAME", RestartFilename, true);

  if( RestartSolution )
    Terminate("CConfig::ReadICOptions", __FILE__, __LINE__,
              "Read restart solution is not supported (yet).");

  // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


bool CConfig::ReadTemporalOptions
(
 const char *configFile
)
 /*
	* Function that reads the temporal specifications.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

	// Read start time of the simulation.
	AddScalarOption(paramFile, "START_TIME", SimulationTime[0], true);
	// Read end time of the simulation.
	AddScalarOption(paramFile, "FINAL_TIME", SimulationTime[1], true);
	// Read time step selected.
	AddScalarOption(paramFile, "TIME_STEP", TimeStep, DefaultParam.TimeStep, true);
	// Read maximum time iteration.
	AddScalarOption(paramFile, "MAX_ITER", MaxIter, DefaultParam.MaxIter, true);
	// Read temporal scheme used in time marching.
	AddScalarOption(paramFile, "TIME_MARCHING", NameTemporalScheme, true);
  // Read input adaptive time step.
  AddBoolOption(paramFile, "ADAPT_TIME", AdaptTime, DefaultParam.NameAdaptTime, true);
  // Read Courant-Friedrichs-Lewy condition number.
  AddScalarOption(paramFile, "CFL_NUMBER", CFL, DefaultParam.CFL, true);
	// Assign TypeTemporalScheme map.
	MapTemporalScheme();

  // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


bool CConfig::ReadBoundaryOptions
(
 const char *configFile
)
 /*
	* Function that reads the boundary information.
	*/
{
  // Open input file.
  std::ifstream paramFile(configFile);

  // Read requested boundary marker, if present.
  AddVectorOption(paramFile, "MARKER_BOUNDARY", NameBoundaryCondition, true);
  // Assign TypeBC map. Note, indices: (SOUTH, NORTH, WEST, EAST).
  MapTypeExternalBC();

  // Read the sponge-layer damping constant.
  AddVectorOption(paramFile, "SPONGE_DAMPING_CONSTANT", DampingConstant, DefaultParam.DampingConstant, true);
  // Read the sponge-damping exponential.
  AddVectorOption(paramFile, "SPONGE_DAMPING_EXPONENT", DampingExponent, DefaultParam.DampingExponent, true);

  // Read the grid-stretching constant.
  AddVectorOption(paramFile, "GRID_STRETCHING_CONSTANT", GridStretchingConstant, DefaultParam.GridStretchingConstant, true);
  // Read the grid-stretching exponential.
  AddVectorOption(paramFile, "GRID_STRETCHING_EXPONENT", GridStretchingExponent, DefaultParam.GridStretchingExponent, true);

  // Close file.
  paramFile.close();

  // Return happily.
	return true;
}


void CConfig::MapTypeExternalBC
(
 void
)
 /*
	* Function that maps NameBoundaryCondition to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["PERIODIC"] = BC_INTERFACE;

  // Initialize actual mapped data.
  TypeExternalBC.resize(nFace);

  for(int iFace=0; iFace<nFace; iFace++){

    // Initialize to unknown.
		TypeExternalBC[iFace] = BC_UNKNOWN;

    // Check if data abides by map convention.
		try {
			Mapper.at(NameBoundaryCondition[iFace]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeExternalBC", __FILE__, __LINE__,
								"BC data does not follow associated map convention!");
		}
    // Assign data according to dedicated enum.
		TypeExternalBC[iFace] = Mapper.at(NameBoundaryCondition[iFace]);
  }
}


void CConfig::MapTypeIC
(
 void
)
 /*
	* Function that maps NameInitialCondition to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["GAUSSIAN_PRESSURE"]   = IC_GAUSSIAN_PRESSURE;
  Mapper["ISENTROPIC_VORTEX"]   = IC_ISENTROPIC_VORTEX;
  Mapper["ENTROPY_WAVE"]        = IC_ENTROPY_WAVE;
  Mapper["VORTEX_ROLLUP"]       = IC_VORTEX_ROLLUP;
  Mapper["ACOUSTIC_PLANE_WAVE"] = IC_ACOUSTIC_PLANE_WAVE;

	// Initialize to unknown.
	TypeIC = IC_UNKNOWN;

	// Check if data abides by map convention.
	try {
		Mapper.at(NameInitialCondition);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTypeIC", __FILE__, __LINE__,
							"IC data does not follow associated map convention!");
	}
	// Assign data according to dedicated enum.
	TypeIC = Mapper.at(NameInitialCondition);
}


void CConfig::MapTemporalScheme
(
 void
)
 /*
	* Function that maps NameTemporalScheme to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["LSRK4"]  = TEMPORAL_SCHEME_LSRK4;
	Mapper["CRK4"]   = TEMPORAL_SCHEME_CRK4;
  Mapper["SSPRK3"] = TEMPORAL_SCHEME_SSPRK3;

	// Initialize to unknown.
	TypeTemporalScheme = TEMPORAL_SCHEME_UNKNOWN;

	// Check if data abides by map convention.
	try {
		Mapper.at(NameTemporalScheme);
	}
	catch ( std::out_of_range& ) {
		Terminate("CConfig::MapTemporalScheme", __FILE__, __LINE__,
							"Temporal data does not follow associated map convention!");
	}

	// Assign data according to dedicated enum.
	TypeTemporalScheme = Mapper.at(NameTemporalScheme);
}


void CConfig::MapTypeStencil
(
 void
)
 /*
	* Function that maps TypeStencil from string to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["DRP_CENTRAL_M3N3"] = STENCIL_DRP_M3N3;
	Mapper["DRP_UPWIND_M2N4"]  = STENCIL_DRP_M2N4;

  // Initialize actual mapped data.
  TypeStencil.resize(nDim);

  // Loop over each dimension.
  for(unsigned short iDim=0; iDim<nDim; iDim++){
  	// Initialize to unknown.
  	TypeStencil[iDim] = STENCIL_UNKNOWN;

  	// Check if data abides by map convention.
  	try {
  		Mapper.at(NameTypeStencil[iDim]);
  	}
  	catch ( std::out_of_range& ) {
  		Terminate("CConfig::MapTypeStencil", __FILE__, __LINE__,
  							"Stencil data does not follow associated map convention!");
  	}

  	// Assign data according to dedicated enum.
  	TypeStencil[iDim] = Mapper.at(NameTypeStencil[iDim]);
  }
}


void CConfig::MapTypeBufferLayer
(
 void
)
 /*
	* Function that maps NameTypeBufferLayer to its enum type.
	*/
{
	// Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

  // Assign mapper according to its enum convention.
  Mapper["NONE"]   = NO_LAYER;
  Mapper["PML"]    = PML_LAYER;
  Mapper["SPONGE"] = SPONGE_LAYER;

  // Number of layers.
  unsigned short nBuffer = NameTypeBufferLayer.size();

  // Initialize to no layer..
	TypeBufferLayer.resize(nBuffer, NO_LAYER);

  // Iterate of data.
	for(int iBuffer=0; iBuffer<nBuffer; iBuffer++){

    // Initialize to no layer.
    TypeBufferLayer[iBuffer] = NO_LAYER;

    // Check if data abides by map convention.
  	try {
  		Mapper.at(NameTypeBufferLayer[iBuffer]);
  	}
  	catch ( std::out_of_range& ) {
  		Terminate("CConfig::MapTypeBufferLayer", __FILE__, __LINE__,
  							"Buffer layer type does not follow associated map convention!");
  	}
    // Assign data according to dedicated enum.
  	TypeBufferLayer[iBuffer] = Mapper.at(NameTypeBufferLayer[iBuffer]);
  }
}


void CConfig::MapTypeZone
(
 void
)
 /*
	* Function that maps TypeZone to its enum type.
	*/
{
	// Initialize actual mapped data.
	TypeZone.resize(nZone);

  // Initialize mapper.
	std::map<std::string, unsigned short> Mapper;

	// Assign mapper according to its enum convention.
	Mapper["ZONE_MAIN"]     = ZONE_MAIN;
	Mapper["ZONE_WEST"]     = ZONE_WEST;
  Mapper["ZONE_EAST"]     = ZONE_EAST;
  Mapper["ZONE_SOUTH"]    = ZONE_SOUTH;
  Mapper["ZONE_NORTH"]    = ZONE_NORTH;
  Mapper["ZONE_CORNER_0"] = ZONE_CORNER_0;
  Mapper["ZONE_CORNER_1"] = ZONE_CORNER_1;
  Mapper["ZONE_CORNER_2"] = ZONE_CORNER_2;
  Mapper["ZONE_CORNER_3"] = ZONE_CORNER_3;

  // Base zone name.
  std::string basename = "ZONE_";

	// Iterate of data.
	for(int iZone=0; iZone<nZone; iZone++){

    // Initialize to unknown.
		TypeZone[iZone] = ZONE_UNKNOWN;

		// Check if data abides by map convention.
		try {
			Mapper.at(NameZoneMarker[iZone]);
		}
		catch ( std::out_of_range& ) {
			Terminate("CConfig::MapTypeZone", __FILE__, __LINE__,
								"Zone names do not follow associated map convention!");
		}
		// Assign data according to dedicated enum.
		TypeZone[iZone] = Mapper.at(NameZoneMarker[iZone]);
	}

	// Make sure input zone markers are unique.
	for(unsigned short jZone=0; jZone<nZone; jZone++)
		for(unsigned short iZone=jZone+1; iZone<nZone; iZone++)
			if( TypeZone[jZone] == TypeZone[iZone] )
				Terminate("CConfig::MapTypeZone", __FILE__, __LINE__,
									"Zone markers must be unique!");
}


void CConfig::DetermineMultizoneStrategy
(
  void
)
 /*
  * Function that determines what strategy for a multizone simulation to adopt.
  */
{
  // Consistency check, input zone markers must be either 1, 2 or 9.
  if( (nZone != 1) && (nZone != 2) && (nZone != 3) && (nZone != 9) )
    Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
              "Multizone strategy demands nZone be: 1, 2, 3 or 9.");

  // To simplify things, make sure the first zone is always the main zone.
  if( TypeZone[0] != ZONE_MAIN )
    Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
              "First input zone must be ZONE_MAIN.");

  // Determine if this is strategy main.
  if( nZone == 1 ){
    // Assign multizone strategy.
    MultizoneStrategy = MULTIZONE_STRATEGY_MAIN;
    return;
  }

  // Determine if this is an full-zonal strategy.
  if( nZone == 9 ){
    // Assign multizone strategy.
    MultizoneStrategy = MULTIZONE_STRATEGY_ALL;
    // Make sure the element ratio input is correct.
    if( (dsNodalRatioZone[0] <= 0.0) || (dsNodalRatioZone.size() != 4) )
      Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                "ELEMENT_RATIO must be positive and: (r(w), r(e), r(s), r(n).");
    return;
   }

   // Determine if this is either combination strategies. Note, the
   // first zone is always fixed as the main one so no need to check for that.
   if( nZone == 2 ){
     switch( TypeZone[1] ){
       case(ZONE_WEST):  MultizoneStrategy = MULTIZONE_STRATEGY_WEST;  break;
       case(ZONE_EAST):  MultizoneStrategy = MULTIZONE_STRATEGY_EAST;  break;
       case(ZONE_SOUTH): MultizoneStrategy = MULTIZONE_STRATEGY_SOUTH; break;
       case(ZONE_NORTH): MultizoneStrategy = MULTIZONE_STRATEGY_NORTH; break;

       default:
         Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                   "Wrong second zone input, must be: WEST, EAST, SOUTH or NORTH.");
     }

     return;
   }

   // Determine if this is either combination strategies. Note, the first zone
   // is always fixed as the main one so no need to check for that.
   if( nZone == 3 ){

     // Make sure all zones are unique.
     if( TypeZone[1] == TypeZone[2] )
      Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                "Second and third zone inputs must be unique.");

     // Check type of second zone.
     switch( TypeZone[1] ){

       // Check horizonal layers.
       case(ZONE_WEST): case(ZONE_EAST):
       {
         if( (TypeZone[2] == ZONE_WEST) || (TypeZone[2] == ZONE_EAST) )
           MultizoneStrategy = MULTIZONE_STRATEGY_HORIZONAL;

         break;
       }

       // Check vertical layers.
       case(ZONE_SOUTH): case(ZONE_NORTH):
       {
         if( (TypeZone[2] == ZONE_SOUTH) || (TypeZone[2] == ZONE_NORTH) )
           MultizoneStrategy = MULTIZONE_STRATEGY_VERTICAL;

         break;
       }

       default:
         Terminate("CConfig::DetermineMultizoneStrategy", __FILE__, __LINE__,
                   "Wrong second/third zone input combinations, must be vertical or horizontal.");
     }

     return;
   }
}


void CConfig::ProcessZoneConformity
(
  void
)
 /*
  * Function that enforces the zone conformity if turned on.
  */
{
  // In case this is a single zone, do nothing and return.
  if( nZone == 1 ) return;

  // Determine what multizone strategy is employed.
  switch(MultizoneStrategy){

    // Combination of east or west zones.
    case(MULTIZONE_STRATEGY_EAST): case(MULTIZONE_STRATEGY_WEST):
    nyNode[1] = nyNode[0]; break;

    // Combination of south or north zones.
    case(MULTIZONE_STRATEGY_SOUTH): case(MULTIZONE_STRATEGY_NORTH):
    nxNode[1] = nxNode[0]; break;

    // All zones are utilized.
    case(MULTIZONE_STRATEGY_ALL):
    {
      // Consistency in naming check.
      assert( TypeZone[0] == ZONE_MAIN );

      for(unsigned short iZone=0; iZone<nZone; iZone++){
        switch(TypeZone[iZone]){
          case(ZONE_MAIN):  break;
          case(ZONE_WEST):  case(ZONE_EAST):  nyNode[iZone] = nyNode[ZONE_MAIN]; break;
          case(ZONE_SOUTH): case(ZONE_NORTH): nxNode[iZone] = nxNode[ZONE_MAIN]; break;
          case(ZONE_CORNER_0):
          {
            nxNode[iZone] = nxNode[ZONE_WEST];
            nyNode[iZone] = nyNode[ZONE_SOUTH];
            break;
          }
          case(ZONE_CORNER_1):
          {
            nxNode[iZone] = nxNode[ZONE_EAST];
            nyNode[iZone] = nyNode[ZONE_SOUTH];
            break;
          }
          case(ZONE_CORNER_2):
          {
            nxNode[iZone] = nxNode[ZONE_WEST];
            nyNode[iZone] = nyNode[ZONE_NORTH];
            break;
          }
          case(ZONE_CORNER_3):
          {
            nxNode[iZone] = nxNode[ZONE_EAST];
            nyNode[iZone] = nyNode[ZONE_NORTH];
            break;
          }
        }
      }

      break;
    }

    // Horizontal zones only.
    case(MULTIZONE_STRATEGY_HORIZONAL):
    {
      nyNode[1] = nyNode[0];
      nyNode[2] = nyNode[0];
      break;
    }

    // Vertical zones only.
    case(MULTIZONE_STRATEGY_VERTICAL):
    {
      nxNode[1] = nxNode[0];
      nxNode[2] = nxNode[0];
      break;
    }

  }
}


void CConfig::dsNodalRatioZone
(
  void
)
 /*
  * Function that checks for inconsistency errors in the nodal ratio
  * specified, w.r.t. the multizone strategy employed.
  */
{
  // In case this is a single-zone simulation, no need to check.
  if( MultizoneStrategy == MULTIZONE_STRATEGY_MAIN ) return;

  // Error flag.
  bool ErrorDetected = false;

  // If this is not a single-zone simulation, make sure all is consistent.
  if( dsNodalRatioZone.size() != nFace ) ErrorDetected = true;
  for(unsigned short i=0; i<nFace; i++)
    if( dsNodalRatioZone[i] <= 0.0)
      ErrorDetected = true;

  // Report error and exit, in case detected.
  if( ErrorDetected )
    Terminate("CConfig::dsNodalRatioZone", __FILE__, __LINE__,
              "NODAL_RATIO is inconsistent: must be positive and of size: 4");
}


