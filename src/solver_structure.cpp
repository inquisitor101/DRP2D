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
  // Number of boundaries in this zone. Note, this is always fixed as 4.
	nBoundary = nFace;
  // Total number of nodes.
  nNode     = geometry_container->GetnNode();
}


CSolver::~CSolver
(
 void
)
 /*
	* Destructor for CSolver class, frees allocated memory.
	*/
{

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
  DataFlux.resize(nDim);
  DataFlux[0].resize(nVar, nullptr);
  DataFlux[1].resize(nVar, nullptr);
  DataSolution.resize(nVar, nullptr);
  DataResidual.resize(nVar, nullptr);

  // Allocate in every variable.
  for(unsigned short iVar=0; iVar<nVar; iVar++){

    // Allocate actual memory.
    DataFlux[0][iVar]  = new as3double[nNode]();
    DataFlux[1][iVar]  = new as3double[nNode]();
    DataSolution[iVar] = new as3double[nNode]();
    DataResidual[iVar] = new as3double[nNode]();

    // Check if allocation failed.
    if( !DataSolution[iVar] || !DataResidual[iVar] || !DataFlux[0][iVar] || !DataFlux[1][iVar] )
      Terminate("CEESolver::CEESolver", __FILE__, __LINE__,
                "Allocation failed for DataSolution, DataResidual or DataFlux.");
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

  for(unsigned short i=0; i<DataFlux.size(); i++)
    if( !DataFlux[i].empty() )
      for(unsigned short j=0; j<DataFlux[i].size(); j++)
        if( DataFlux[i][j] ) delete [] DataFlux[i][j];
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
  }
}


as3double CEESolver::ComputeTimeStep
(
 void
)
 /*
  * Function that computes the time step per input element.
  */
{
  // Time step -- dummy for now.
  as3double dt = 0.0;

  // Return computed value.
  return dt;
}