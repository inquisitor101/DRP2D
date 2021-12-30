#include "driver_structure.hpp"



CDriver::CDriver
(
 CConfig *config
)
 /*
	* Constructor, used to set-up, run and (post-)process data.
	*/
{
	// Initialize all containers to nullptr.
	SetContainers_Null();

  // Choose active config container as input container.
  config_container = config;

  // Preprocess geometry container.
  Geometry_Preprocessing(config_container);

  // Preprocess input container.
  Input_Preprocessing(config_container,
  										geometry_container);

  // Preprocess output container.
  Output_Preprocessing(config_container,
  										 geometry_container);

  // Preprocess initial container.
  Initial_Preprocessing(config_container,
  											geometry_container);

  // Preprocess spatial container.
  Spatial_Preprocessing(config_container,
  											geometry_container,
                        initial_container);

  // Preprocess solver container.
  Solver_Preprocessing(config_container,
  										 geometry_container,
                       initial_container,
  										 spatial_container);

  // Preprocess iteration container.
  Iteration_Preprocessing(config_container,
  												geometry_container,
  												solver_container,
  												spatial_container);

  // Preprocess the parallelization. Note, this must be set before the
  // temporal container.
  Parallelization_Preprocessing();

  // Preprocess temporal container.
  Temporal_Preprocessing(config_container,
  											 geometry_container,
  											 iteration_container,
  											 solver_container,
  											 spatial_container);


  // Simulation start time.
  SimTimeStart = config_container->GetSimulationStartTime();
  // Simulation end time.
  SimTimeFinal = config_container->GetSimulationFinalTime();

  // Max temporal iterations specified.
  MaxTimeIter  = config_container->GetMaxIter();
}


CDriver::~CDriver
(
 void
)
 /*
	* Destructor for CDriver class, frees allocated memory.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Freeing memory..." << std::endl;

	std::cout << "  deleting CConfig........ ";
	if(config_container != nullptr) delete config_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CGeometry...... ";
	if(geometry_container != nullptr) delete geometry_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CTemporal...... ";
	if(temporal_container != nullptr) delete temporal_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CInput......... ";
	if(input_container != nullptr) delete input_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting COutput........ ";
	if(output_container != nullptr) delete output_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CIteration..... ";
	if(iteration_container != nullptr) delete iteration_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CSolver........ ";
	if(solver_container != nullptr) delete solver_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CInitial....... ";
	if(initial_container != nullptr) delete initial_container;
	std::cout << "Done." << std::endl;

	std::cout << "  deleting CSpatial....... ";
	if(spatial_container != nullptr) delete spatial_container;
	std::cout << "Done." << std::endl;

  std::cout << "  deleting CProcess....... ";
	if(process_container != nullptr) delete process_container;
	std::cout << "Done." << std::endl;

	// Set all containers deleted to nullptr.
	config_container    = nullptr;
	geometry_container  = nullptr;
	input_container     = nullptr;
	output_container    = nullptr;
	solver_container    = nullptr;
	temporal_container  = nullptr;
	iteration_container = nullptr;
	spatial_container   = nullptr;
	initial_container   = nullptr;
  process_container   = nullptr;

	std::cout << "Done." << std::endl;
}


void CDriver::Parallelization_Preprocessing
(
 void
)
 /*
  * Function that sets the necessary steps needed for a parallel implementation.
  */
{
	// Total number of nodes.
  unsigned long nNode = geometry_container->GetnNode();

  // Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
#ifdef HAVE_OPENMP
  // Get max number of threads specified.
  unsigned short nThreads = omp_get_max_threads();
  std::cout << "This is a parallel implementation using: "
            << nThreads << " threads." << std::endl;

  // Check if the work load is divisible by the number of threads.
  if( nNode%nThreads != 0 )
    std::cout << "\n**********************************************"
  						<< "**********************************************\n"
              << "Warning, inefficient implementation!\n"
              << "... Number of total nodes is not a multiple of the number of threads."
              << "\n**********************************************"
              << "**********************************************\n" << std::endl;

  // Estimate computational work load of each thread in terms of nodes.
  unsigned long WorkLoadNode = nNode/omp_get_max_threads();

  // Report work load per element.
  std::cout << "Workload of nodes shared among each thread is: "
            << WorkLoadNode << " [Node/Thread]. " << std::endl;

#else
  std::cout << "This is a serial implementation." << std::endl;
#endif

}


void CDriver::SetContainers_Null
(
 void
)
 /*
	* Function that initializes all containers to nullptr.
	*/
{
	// Containers initialized.
	config_container    = nullptr;
	geometry_container  = nullptr;
	temporal_container  = nullptr;
	input_container     = nullptr;
	output_container    = nullptr;
	solver_container    = nullptr;
	iteration_container = nullptr;
	spatial_container   = nullptr;
	initial_container   = nullptr;
  process_container   = nullptr;
}


void CDriver::Geometry_Preprocessing
(
 CConfig *config_container
)
 /*
	* Function that preprocesses the geometry container.
	*/
{
	// Assign a geometry container.
	geometry_container = new CGeometry(config_container);
}


void CDriver::Input_Preprocessing
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that preprocesses the input container.
	*/
{
	// Assign an input container.
  input_container = new CInput(config_container,
									 						 geometry_container);
}


void CDriver::Output_Preprocessing
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that preprocesses the output container.
	*/
{
	// Assign an output container.
	output_container = new COutput(config_container,
																 geometry_container);
}


void CDriver::Iteration_Preprocessing
(
 CConfig     *config_container,
 CGeometry   *geometry_container,
 CSolver     *solver_container,
 CSpatial    *spatial_container
)
 /*
	* Function that preprocesses the iteration container.
	*/
{
  // Check which iteration to specify.
  switch( config_container->GetTypeSolver() ){

  	// Euler solver.
  	case(SOLVER_EE):
  	{
      // Assign iteration container.
  		iteration_container = new CEEIteration(config_container,
  																				 	 geometry_container,
  																					 solver_container);
  		break;
  	}

  	// Otherwise, exit immediately.
  	default:
  		Terminate("CDriver::Iteration_Preprocessing", __FILE__, __LINE__,
  							"Iteration strategy for solver specified is not (yet) implemented!");
  }
}


void CDriver::Solver_Preprocessing
(
 CConfig   	 *config_container,
 CGeometry 	 *geometry_container,
 CInitial    *initial_container,
 CSpatial    *spatial_container
)
 /*
	* Function that preprocesses the solver container.
	*/
{
  // Check which solver to specify.
  switch( config_container->GetTypeSolver() ){

  	// Euler solver.
  	case(SOLVER_EE):
  	{
      // Assign solver container.
  		solver_container = new CEESolver(config_container,
																			 geometry_container,
                                       initial_container,
																			 spatial_container);
  		break;
  	}

  	// Otherwise, exit immediately.
  	default:
  		Terminate("CDriver::Solver_Preprocessing", __FILE__, __LINE__,
  							"Solver container specified is not (yet) implemented!");
  }
}


void CDriver::Spatial_Preprocessing
(
 CConfig   	*config_container,
 CGeometry 	*geometry_container,
 CInitial   *initial_container
)
 /*
	* Function that preprocesses the spatial discretization container.
	*/
{
  // Check which spatial discretization to specify.
  switch( config_container->GetTypeSolver() ){

    // Euler solver.
    case(SOLVER_EE):
    {
      // Assign spatial class.
      spatial_container = new CEESpatial(config_container,
    																	   geometry_container,
                                         initial_container);
      break;
    }

    // Otherwise, exit immediately.
    default:
      Terminate("CDriver::Spatial_Preprocessing", __FILE__, __LINE__,
                "Spatial container for solver specified is not (yet) implemented!");

	}
}


void CDriver::Process_Preprocessing
(
 CConfig    *config_container,
 CGeometry  *geometry_container,
 COutput    *output_container,
 CInitial   *initial_container,
 CSolver    *solver_container,
 CSpatial   *spatial_container
)
 /*
  * Function that preprocesses the process container.
  */
{
  // Check which process to specify.
  switch( config_container->GetTypeSolver() ){

    // Pure EE solver.
    case(SOLVER_EE):
    {
      // Assign process container.
      process_container = new CEEProcess(config_container,
                                         geometry_container,
                                         output_container,
                                         initial_container,
                                         solver_container,
                                         spatial_container);
      break;
    }

    // Otherwise, exit immediately.
    default:
      Terminate("CDriver::Process_Preprocessing", __FILE__, __LINE__,
                "Process container for specified solver is not (yet) implemented!");
  }
}


void CDriver::Initial_Preprocessing
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that preprocesses the initial condition container.
	*/
{
  // Check which initial condition to specify.
  switch( config_container->GetTypeIC() ){

  	// Gaussian IC.
  	case(IC_GAUSSIAN_PRESSURE):
  	{
      // Assign initial container.
  		initial_container = new CGaussianInitial(config_container,
																				  		 geometry_container);
  		break;
  	}

    // Isentropic vortex IC.
    case(IC_ISENTROPIC_VORTEX):
    {
      // Assign initial container.
      initial_container = new CIsentropicVortexInitial(config_container,
                                                       geometry_container);
      break;
    }

  	// Otherwise, exit immediately.
  	default:
  		Terminate("CDriver::Initial_Preprocessing", __FILE__, __LINE__,
  							"Initial condition prescribed is not (yet) implemented!");
  }
}


void CDriver::Temporal_Preprocessing
(
 CConfig     *config_container,
 CGeometry   *geometry_container,
 CIteration  *iteration_container,
 CSolver     *solver_container,
 CSpatial    *spatial_container
)
 /*
	* Function that preprocesses the temporal container.
	*/
{
	// Assign temporal container.
	switch( config_container->GetTypeTemporalScheme() ){

		// Low-storage 4th-order 5-stage Runge-Kutta (explicit).
		case(TEMPORAL_SCHEME_LSRK4):
		{
      // Assign temporal container.
			temporal_container = new CLSRK4Temporal(config_container,
																							geometry_container,
																							iteration_container,
																							solver_container,
																							spatial_container);
			break;
		}

    // // Strong-stability-preserving 3-stage Runge-Kutta (explicit).
		// case(TEMPORAL_SCHEME_SSPRK3):
		// {
    //   // Assign temporal container.
		// 	temporal_container = new CSSPRK3Temporal(config_container,
		// 																					 geometry_container,
		// 																					 iteration_container,
		// 																					 solver_container,
		// 																					 spatial_container);
		// 	break;
		// }

    // Otherwise, exit immediately.
		default:
			Terminate("CDriver::Temporal_Preprocessing", __FILE__, __LINE__,
								"Temporal scheme specified is not (yet) implemented!");
	}
}


void CDriver::MonitorOutput
(
 unsigned long           iIter,
 as3double               time,
 as3double               dt,
 as3vector1d<as3double> &MonitoringData,
 bool                    MonitorData = true
)
 /*
	* Function that outputs the header of the information being displayed.
	*/
{
	// Number of output reports for monitoring progress.
	unsigned long nOutput = std::max(1ul, MaxTimeIter/100);
	// Compute number of max digits needed for output.
	unsigned long nDigits = std::to_string(MaxTimeIter).size();
  // Monitoring output frequency.
  unsigned long OutputFreq = config_container->GetOutputFreq();

	// Display header.
	if( iIter%nOutput == 0 ){
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
		std::cout << " Iteration\tPhysical Time \t Time step \t Max(Mach) \t Res[RMS(rho)]" << std::endl;
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
	}

	// Display data in this iteration.
	if( MonitorData && (iIter%OutputFreq==0) ){

    // Extract the maximum Mach number.
    const as3double Mmax = MonitoringData[0];

		// Display progress.
		std::cout << std::scientific << "   "
							<< std::setw(nDigits) << iIter
							<< " \t "  << time
							<< " \t "  << dt
              << " \t "  << Mmax
							<< " \t "  << "N/A" << std::endl;
	}


  // Check for floating-point errors at run-time.
#ifdef ENABLE_NAN_CHECK
    CheckFloatingError();
#endif
}


void CDriver::StartSolver
(
 void
)
 /*
	* Function that initiates the solver driver.
	*/
{
	// Gauge start time.
	// as3double startTime = as3double(clock())/as3double(CLOCKS_PER_SEC);
  as3double startTime = omp_get_wtime();

	// Run a preprocessing step to initialize the solution and condition the data
	// in case there need be. Note, this in only executed once and before marching
	// in time.
	Preprocess();

	// Begin actual solver.
	Run();

	// Gauge end time used by the solver.
	// as3double stopTime = as3double(clock())/as3double(CLOCKS_PER_SEC);
  as3double stopTime = omp_get_wtime();

	// Lapse time used by the entire solver.
	as3double lapsedTime = stopTime - startTime;


	// Report lapsed time.
	std::cout << "\n% % % % % % % % % % % % % % % % %" << std::endl;
	std::cout << std::scientific << "lapsed time [sec]: " << lapsedTime << " %" << std::endl;
}


void CDriver::Preprocess
(
 void
)
 /*
	* Function that preprocesses the solver, before commencing with the iterative
	* solution.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Initializing solution... ";


  // Initialize solution.
  solver_container->InitializeSolution(geometry_container,
                                       initial_container,
                                       0.0);

  // Report progress.
	std::cout << "Done." << std::endl;

	// Write output VTK for initial condition.
	output_container->WriteFileVTK(config_container,
																 geometry_container,
																 solver_container);
}


void CDriver::Run
(
 void
)
 /*
	* Function that runs the entire solver.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
	switch( config_container->GetTypeSolver() ){
		case(SOLVER_EE): std::cout << "Beginning EE Solver." << std::endl; break;
	}
	std::cout << std::endl;

  // Output writing frequency.
  const unsigned long WriteFreq = config_container->GetWriteFreq();

  // Monitoring data.
  // Thus far, use only [0]: max(Mach).
  as3vector1d<as3double> MonitoringData(1);

	// Estimate time step needed.
	as3double dt = solver_container->ComputeTimeStep(config_container, geometry_container);

	// Current time.
	as3double SimTime = SimTimeStart;
	// Current iteration.
	unsigned long IterCount = 0;
  // Marker to check if final time step is written or not.
  bool FinalStep = false;


	// Display header for output format.
	MonitorOutput(IterCount, SimTime, dt, MonitoringData, false);

	// March in time, until target time is reached.
	while( (SimTime < SimTimeFinal) && (IterCount < MaxTimeIter) ){

		// Execute a single temporal update.
		temporal_container->TimeMarch(config_container,
																	geometry_container,
																	iteration_container,
																	solver_container,
																	spatial_container,
                                  initial_container,
																	SimTime, dt,
                                  MonitoringData);

		// Update (physical) time.
		SimTime += dt;

		// Update iteration count.
		IterCount++;


    // Check if this is the final time step.
    if( (SimTime >= SimTimeFinal) || (IterCount >= MaxTimeIter) )
      FinalStep = true;

    // Process data every OutputFreq iterations.
    if( IterCount%WriteFreq == 0 || FinalStep ){

      // Write output VTK for initial condition.
      output_container->WriteFileVTK(config_container,
    																 geometry_container,
    																 solver_container);

      // If this is an adaptive time-stepping, then compute new time-step.
      if( config_container->GetAdaptTime() )
        dt = solver_container->ComputeTimeStep(config_container, geometry_container);
    }

		// Display output for progress monitoring.
		MonitorOutput(IterCount, SimTime, dt, MonitoringData);
	}
}