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

  // Consistency check.
  if( NameTypeStencil.size() != nDim )
    Terminate("CConfig::ReadSolverOptions", __FILE__, __LINE__,
              "TYPE_STENCIL must be of dimension: 2.");

  // Deduce stencil dimensions.
  for(unsigned short iDim=0; iDim<nDim; iDim++){
    switch( TypeStencil[iDim] ){
      case(STENCIL_DRP_M3N3): NStencil[iDim] = 3; MStencil[iDim] = 3; break;
      case(STENCIL_DRP_M2N4): NStencil[iDim] = 4; MStencil[iDim] = 2; break;
      default:
        Terminate("CConfig::ReadSolverOptions", __FILE__, __LINE__,
                  "Type of stencil is not supported (yet).");
    }
  }

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

	// Read number of nodes in x-direction.
	AddScalarOption(paramFile, "NUMBER_XNODE", nxNode, true);
	// Read number of nodes in y-direction.
	AddScalarOption(paramFile, "NUMBER_YNODE", nyNode, true);

  // Read domain bounds in each zone.
	AddVectorOption(paramFile, "DOMAIN_BOUND", DomainBound, true);
  // Read buffer layer starting location in x-direction.
  AddScalarOption(paramFile, "BUFFER_INTERFACE_X", xBufferInterface, true);

  // Assign solver type. For now, use a EE only.
  TypeSolver = SOLVER_EE;

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

  // Read damping constant.
  AddScalarOption(paramFile, "SPONGE_DAMPING_CONSTANT", DampingConstant, DefaultParam.DampingConstant, true);
  // Read damping exponential.
  AddScalarOption(paramFile, "SPONGE_DAMPING_EXPONENT", DampingExponent, DefaultParam.DampingExponent, true);

  // Read the grid-stretching constant.
  AddScalarOption(paramFile, "GRID_STRETCHING_CONSTANT", GridStretchingConstant, DefaultParam.GridStretchingConstant, true);
  // Read the grid-stretching exponential.
  AddScalarOption(paramFile, "GRID_STRETCHING_EXPONENT", GridStretchingExponent, DefaultParam.GridStretchingExponent, true);

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




