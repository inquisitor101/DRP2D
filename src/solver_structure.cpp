#include "solver_structure.hpp"




CSolver::CSolver
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CInitial  *initial_container,
 CSpatial  *spatial_container
)
 /*
	* Constructor, used to initialize CSolver in zone: iZone.
	*/
{
  // Total number of nodes.
  nNode     = geometry_container->GetnNode();
  // Number of boundaries in this zone. Note, this is always fixed as 4.
	nBoundary = nFace;

  // Initialize boundary data container.
  Boundary_Preprocessing(config_container,
                         geometry_container,
                         initial_container);
}


CSolver::~CSolver
(
 void
)
 /*
	* Destructor for CSolver class, frees allocated memory.
	*/
{
  for(unsigned short i=0; i<boundary_container.size(); i++)
  	if( boundary_container[i] ) delete boundary_container[i];
}


void CSolver::Boundary_Preprocessing
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CInitial  *initial_container
)
 /*
	* Function that preprocess the boundary container and initialize it.
	*/
{
  // Header for output.
  std::cout << "----------------------------------------------"
               "----------------------------------------------\n";
  std::cout << "Processing boundary of type: "
            << DisplaySolverType( config_container->GetTypeSolver() )
            << std::endl;

	// Initialize the needed number of boundaries in this zone.
	boundary_container.resize(nBoundary, nullptr);

	// Loop over every boundary and specify input condition.
	for(unsigned short iBoundary=0; iBoundary<nBoundary; iBoundary++){

		// Check which boundary condition to apply.
		switch( config_container->GetTypeExternalBC(iBoundary) ){

			case(BC_INTERFACE):
			{
				// Initialize interface/periodic boundary.
        switch( config_container->GetTypeSolver() ){

          // Pure EE interface boundary.
          case( SOLVER_EE ):
          {

            // Initialize interface/periodic boundary.
            boundary_container[iBoundary] = new CEEInterfaceBoundary(config_container,
                                                                     geometry_container,
                                                                     iBoundary);
            break;
          }

          default:
            Terminate("CSolver::Boundary_Preprocessing", __FILE__, __LINE__,
                      "Unknown solver specified for the interface boundary.");
        }

        // Break out of the interface case.
        break;
			}

			default:
				Terminate("CSolver::Boundary_Preprocessing", __FILE__, __LINE__,
									"Boundary condition is unknown!");
		}
	}

  // Report progress.
  std::cout << "Done." << std::endl;
}


CEESolver::CEESolver
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CInitial  *initial_container,
 CSpatial  *spatial_container
)
	:
		CSolver
		(
		 config_container,
		 geometry_container,
     initial_container,
		 spatial_container
		)
 /*
	* Constructor, used to initialize CEESolver in zone: iZone.
	*/
{
  // Reserve memory for data.
  DataSolution.resize(nVar, nullptr);
  DataResidual.resize(nVar, nullptr);

  // Allocate in every variable.
  for(unsigned short iVar=0; iVar<nVar; iVar++){

    // Allocate actual memory.
    DataSolution[iVar] = new as3double[nNode]();
    DataResidual[iVar] = new as3double[nNode]();

    // Check if allocation failed.
    if( !DataSolution[iVar] || !DataResidual[iVar] )
      Terminate("CEESolver::CEESolver", __FILE__, __LINE__,
                "Allocation failed for DataSolution or DataResidual.");
  }
}


CEESolver::~CEESolver
(
 void
)
 /*
	* Destructor for CEESolver class, frees allocated memory.
	*/
{
  for(unsigned short i=0; i<DataSolution.size(); i++)
    if( DataSolution[i] ) delete [] DataSolution[i];

  for(unsigned short i=0; i<DataResidual.size(); i++)
    if( DataResidual[i] ) delete [] DataResidual[i];
}


void CEESolver::InitializeSolution
(
 CGeometry *geometry_container,
 CInitial  *initial_container,
 as3double  time
)
 /*
  * Function that initializes the solution.
  */
{
  // Abbreviation.
  const as3double ovgm1 = 1.0/GAMMA_MINUS_ONE;

  // Extract grid coordinates.
  const as3double *xcoord = geometry_container->GetGridCoordinate(0);
  const as3double *ycoord = geometry_container->GetGridCoordinate(1);

  // Temporary storage for nodal IC.
  as3vector1d<as3double> Q(nVar);

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), private(Q)
#endif
  // Loop over all the nodes.
  for(unsigned long l=0; l<nNode; l++){

    // Explicitly extract current coordinate.
    const as3double x = xcoord[l];
    const as3double y = ycoord[l];

    // Extract the current IC in terms of primitive variables.
    Q = initial_container->SetInitialConditionPrimitive(x, y, time);

    // Explicitly extract primitive data.
    const as3double rho = Q[0];
    const as3double u   = Q[1];
    const as3double v   = Q[2];
    const as3double p   = Q[3];

    // Compute the energy.
    const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

    // Assemble working solution in conervative form.
    DataSolution[0][l] = rho;
    DataSolution[1][l] = rho*u;
    DataSolution[2][l] = rho*v;
    DataSolution[3][l] = rhoE;

  } // End of OpenMP parallel region.
}


as3double CEESolver::ComputeTimeStep
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
  * Function that computes the time step per input element.
  */
{
  // Extract prescribed CFL number.
  const as3double CFL = config_container->GetCFL();

  // Check whether or not a fixed time step is specified.
  if( CFL < 0.0 ){

    // Make sure no adaptive time-stepping is used.
    if( config_container->GetAdaptTime() )
      Terminate("CEESolver::ComputeTimeStep", __FILE__, __LINE__,
                "This is a fixed-time-step simulation, there should be no adaptive time-stepping specified.");

    // Specify the time step input.
    const as3double dt = config_container->GetTimeStep();

    // Return data.
    return dt;
  }

  // Abbreviations.
  const as3double gm1 = GAMMA - 1.0;

  // Extract grid resolution.
  const as3double dx = geometry_container->GetGridResolution(0);
  const as3double dy = geometry_container->GetGridResolution(1);

  // Time step and max resolution-per-eigenvalue ratio.
  as3double dt = 0.0, minratio = 1.0e5;

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(min:minratio)
#endif
  // Loop over all nodes and determine max time step.
  for(unsigned long l=0; l<nNode; l++){

    // Abbreviation: inverse of density.
    const as3double ovrho = 1.0/DataSolution[0][l];
    // Extract primitive variables.
    const as3double u     = ovrho*DataSolution[1][l];
    const as3double v     = ovrho*DataSolution[2][l];
    const as3double p     = gm1*(DataSolution[3][l]
                          - 0.5*(u*DataSolution[1][l] + v*DataSolution[2][l]) );

    // Compute local speed-of-sound.
    const as3double a     = sqrt(GAMMA*p*ovrho);
    // Determine max ratio between grid-resolution per eigenvalue.
    const as3double ratio = std::min( fabs(u)/dx, fabs(v)/dy ) + a;

    // Compute maximum eigenvalue.
    minratio = std::min( minratio, ratio );

  } // End of OpenMP parallel loop.

  std::cout << minratio << std::endl;

  // Determine max allowable time step.
  dt = CFL/minratio;

  // Return computed value.
  return dt;
}



