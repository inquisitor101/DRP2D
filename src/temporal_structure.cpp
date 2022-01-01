#include "temporal_structure.hpp"



CTemporal::CTemporal
(
 CConfig     *config_container,
 CGeometry   *geometry_container,
 CIteration  *iteration_container,
 CSolver     *solver_container,
 CSpatial    *spatial_container
)
 /*
	* Constructor, used to initialize CTemporal.
	*/
{
  // Extract number of nodes in x-direction.
  nxNode = config_container->GetnxNode();
  // Extract number of nodes in y-direction.
  nyNode = config_container->GetnyNode();
  // Total number of nodes in grid.
  nNode  = geometry_container->GetnNode();

  // Extract stencil widths, per dimension.
  MxStencil = config_container->GetmxStencil();
  NxStencil = config_container->GetnxStencil();
  MyStencil = config_container->GetmyStencil();
  NyStencil = config_container->GetnyStencil();

  // Total nodes in x-diretcion.
  nxTotal = nxNode + MxStencil + NxStencil;
  // Total nodes in y-diretcion.
  nyTotal = nyNode + MyStencil + NyStencil;

  // Initialize working array parameter dimensions.
  InitializeWorkArrayDimension(config_container);
  // Reserve needed memory for the working array.
  InitializeWorkArray();
}


CTemporal::~CTemporal
(
 void
)
 /*
	* Destructor for CTemporal class, frees allocated memory.
	*/
{
  for(unsigned short i=0; i<work_array.size(); i++)
    if( work_array[i] ) delete [] work_array[i];
}


void CTemporal::InitializeWorkArrayDimension
(
 CConfig *config_container
)
 /*
  * Function that defines the dimension parameters used in the working array.
  */
{
  // Modified (inclusive ghost nodes) nodal dimension in x-direction.
  unsigned long nx  = nxNode + MxStencil + NxStencil;
  // Modified (inclusive ghost nodes) nodal dimension in y-direction.
  unsigned long ny  = nyNode + MyStencil + NyStencil;
  // Number of grid DOFs for each data value, including periodic ghost nodes.
  nWorkingArrayDOFs = nx*ny;
  // Number of variables per each working array entry.
  nWorkingArrayVar  = nVar;

  // Determine number of entries needed for the working array.
  // Entries: 2,
  // [0]: x-flux pre-computation.
  // [1]: y-flux pre-computation.
  nWorkingArrayEntry = 2;
}


void CTemporal::InitializeWorkArray
(
 void
)
 /*
  * Function that reserves and initializes the required memory for a work array.
  */
{
  // Total number of data needed for 1st/outer index in the working array.
  const unsigned short nDataOuter = nWorkingArrayEntry*nWorkingArrayVar;
  // Total number of data needed for 2nd/inner index in the working array.
  const unsigned long  nDataInner = nWorkingArrayDOFs;

  // Initialize the work array needed to carry out a residual update over an
  // entire grid sweep iteration.
  work_array.resize(nDataOuter, nullptr);

  // Allocate data per every entry of the working array.
  for(unsigned short i=0; i<work_array.size(); i++){

    // Allocate the actual memory.
    work_array[i] = new as3double[nDataInner]();

    // Check if allocation failed.
    if( !work_array[i] )
      Terminate("CTemporal::InitializeWorkArray", __FILE__, __LINE__,
                "Allocation failed for work_array.");
  }
}


CLSRK4Temporal::CLSRK4Temporal
(
 CConfig     *config_container,
 CGeometry   *geometry_container,
 CIteration  *iteration_container,
 CSolver     *solver_container,
 CSpatial    *spatial_container
)
	:
		CTemporal
		(
		 config_container,
		 geometry_container,
		 iteration_container,
		 solver_container,
		 spatial_container
		)
 /*
	* Constructor, used to initialize CLSRK4Temporal.
	*/
{
	// Initialize the LSRK4 coefficients.
	rk4a.resize(nStageRK);
	rk4b.resize(nStageRK);
	rk4c.resize(nStageRK);

	// Low-storage 4th-order Runge-Kutta coefficients.
	rk4a[0] =  0.0;
	rk4a[1] = -567301805773.0/1357537059087.0;
	rk4a[2] = -2404267990393.0/2016746695238.0;
	rk4a[3] = -3550918686646.0/2091501179385.0;
	rk4a[4] = -1275806237668.0/842570457699.0;

	rk4b[0] =  1432997174477.0/9575080441755.0;
	rk4b[1] =  5161836677717.0/13612068292357.0;
	rk4b[2] =  1720146321549.0/2090206949498.0;
	rk4b[3] =  3134564353537.0/4481467310338.0;
	rk4b[4] =  2277821191437.0/14882151754819.0;

	rk4c[0] =  0.0;
	rk4c[1] =  1432997174477.0/9575080441755.0;
	rk4c[2] =  2526269341429.0/6820363962896.0;
	rk4c[3] =  2006345519317.0/3224310063776.0;
	rk4c[4] =  2802321613138.0/2924317926251.0;

  // Reserve memory for tentative data.
  DataSolutionTentative.resize(nVar, nullptr);

  // Allocate in every variable.
  for(unsigned short iVar=0; iVar<nVar; iVar++){

    // Allocate actual memory.
    DataSolutionTentative[iVar] = new as3double[nNode]();

    // Check if allocation failed.
    if( !DataSolutionTentative[iVar] )
      Terminate("CLSRK4Temporal::CLSRK4Temporal", __FILE__, __LINE__,
                "Allocation failed for DataSolutionTentative.");
  }
}


CLSRK4Temporal::~CLSRK4Temporal
(
 void
)
 /*
	* Destructor for CLSRK4Temporal class, frees allocated memory.
	*/
{
  for(unsigned short i=0; i<DataSolutionTentative.size(); i++)
    if( DataSolutionTentative[i] ) delete [] DataSolutionTentative[i];
}


void CLSRK4Temporal::TimeMarch
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration             *iteration_container,
 CSolver                *solver_container,
 CSpatial               *spatial_container,
 CInitial               *initial_container,
 as3double               physicalTime,
 as3double               dtTime,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that updates the simulation by a single LSRK4 step.
	*/
{
	// Loop over all RK stages.
	for(unsigned short iStageRK=0; iStageRK<nStageRK; iStageRK++){

		// Local physical time.
		const as3double localTime = physicalTime + rk4c[iStageRK]*dtTime;

		// Local RK coefficients.
		const as3double alpha = rk4a[iStageRK];
		const as3double beta  = rk4b[iStageRK];

    // Perform an entire LSRK4 sweep and update the residual.
    UpdateTime(config_container,
               geometry_container,
               iteration_container,
               solver_container,
               spatial_container,
               initial_container,
               localTime, dtTime,
               alpha, beta,
               MonitoringData);
	}
}


void CLSRK4Temporal::UpdateTime
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration             *iteration_container,
 CSolver                *solver_container,
 CSpatial               *spatial_container,
 CInitial               *initial_container,
 as3double               localTime,
 as3double               dtTime,
 as3double               alpha,
 as3double               beta,
 as3vector1d<as3double> &MonitoringData

)
 /*
	* Function that performs a single LSRK4 stage sweep and update the residual.
	*/
{
  // Initialize max of the Mach number squared.
  as3double M2max = 0.0;

  // Extract current total solution.
  auto& sol = solver_container->GetDataSolution();
  // Extract current total residual.
  auto& res = solver_container->GetDataResidual();
  // Extract tentative solution.
  auto& tmp = DataSolutionTentative;

  // Total solution size in first dimension.
  unsigned short nSolSize = sol.size();

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

    // Loop over all nodes and copy solution to physical nodes of working array.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3)
#endif
    // Loop over every variable and every node.
    for(unsigned short iVar=0; iVar<nSolSize; iVar++){
      for(unsigned long j=0; j<nyNode; j++){
        for(unsigned long i=0; i<nxNode; i++){
          // Ghost node index.
          unsigned long I = i + NxStencil;
          unsigned long J = j + NyStencil;

          // Equivalent 1D running index of ghost node.
          unsigned long IJ = J*nxTotal + I;
          // Equivalent 1D running index of physical node.
          unsigned long ij = j*nxNode  + i;

          // Copy data point.
          work_array[iVar][IJ] = sol[iVar][ij];
        }
      }
    }


    // Impose the boundary conditions on all boundary nodes.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned short iBoundary=0; iBoundary<nFace; iBoundary++){

      // Extract the relevant boundary container.
      auto* boundary_container = solver_container->GetBoundaryContainer(iBoundary);

      // Impose boundary condition.
      boundary_container->ImposeBoundaryCondition(config_container,
                                                  geometry_container,
                                                  work_array,
                                                  localTime);
    }

  } // End of OpenMP parallel region.


  // Preprocess the data, if needed.
  iteration_container->Preprocess(config_container,
                                  geometry_container,
                                  solver_container,
                                  spatial_container,
                                  work_array,
                                  localTime);

  // Compute the residual on all physical nodes.
  spatial_container->ComputeResidual(config_container,
                                     geometry_container,
                                     spatial_container,
                                     initial_container,
                                     work_array,
                                     res,
                                     localTime,
                                     MonitoringData);

  // Set the max of the Mach squared in this element.
  M2max = MonitoringData[0];


  // Loop over all nodes and update the solution.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(2)
#endif
  // Loop over every variable and node.
  for(unsigned short iVar=0; iVar<nSolSize; iVar++){
    for(unsigned long l=0; l<nNode; l++){
      // Predictor step, update the tentative solution.
      tmp[iVar][l]  = alpha*tmp[iVar][l] + dtTime*res[iVar][l];
      // Corrector step, update the true solution.
      sol[iVar][l] += beta*tmp[iVar][l];
    }
  } // End of OpenMP parallel region.


  // Assign the actual max of the Mach number.
  MonitoringData[0] = sqrt(M2max);
}


CSSPRK3Temporal::CSSPRK3Temporal
(
 CConfig     *config_container,
 CGeometry   *geometry_container,
 CIteration  *iteration_container,
 CSolver     *solver_container,
 CSpatial    *spatial_container
)
	:
		CTemporal
		(
		 config_container,
		 geometry_container,
		 iteration_container,
		 solver_container,
		 spatial_container
		)
 /*
	* Constructor, used to initialize CSSPRK3Temporal.
	*/
{
	// Initialize the LSRK4 coefficients.
	rk4a.resize(nStageRK);
	rk4b.resize(nStageRK);
	rk4c.resize(nStageRK);

	// Strong-stability-preserving 3rd-order Runge-Kutta coefficients.
	rk4a[0] = 1.0;
	rk4a[1] = 3.0/4.0;
	rk4a[2] = 1.0/3.0;

	rk4b[0] = 1.0;
	rk4b[1] = 1.0/4.0;
	rk4b[2] = 2.0/3.0;

	rk4c[0] = 0.0;
	rk4c[1] = 1.0;
	rk4c[2] = 1.0/2.0;

  // Reserve memory for tentative data.
  DataSolutionTentative.resize(nVar, nullptr);

  // Allocate in every variable.
  for(unsigned short iVar=0; iVar<nVar; iVar++){

    // Allocate actual memory.
    DataSolutionTentative[iVar] = new as3double[nNode]();

    // Check if allocation failed.
    if( !DataSolutionTentative[iVar] )
      Terminate("CSSPRK3Temporal::CSSPRK3Temporal", __FILE__, __LINE__,
                "Allocation failed for DataSolutionTentative.");
  }
}


CSSPRK3Temporal::~CSSPRK3Temporal
(
 void
)
 /*
	* Destructor for CSSPRK3Temporal class, frees allocated memory.
	*/
{
  for(unsigned short i=0; i<DataSolutionTentative.size(); i++)
    if( DataSolutionTentative[i] ) delete [] DataSolutionTentative[i];
}


void CSSPRK3Temporal::TimeMarch
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration             *iteration_container,
 CSolver                *solver_container,
 CSpatial               *spatial_container,
 CInitial               *initial_container,
 as3double               physicalTime,
 as3double               dtTime,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that updates the simulation by a single SSPRK3 step.
	*/
{
	// Loop over all RK stages.
	for(unsigned short iStageRK=0; iStageRK<nStageRK; iStageRK++){

		// Local physical time.
		const as3double localTime = physicalTime + rk4c[iStageRK]*dtTime;

		// Local RK coefficients.
		const as3double alpha = rk4a[iStageRK];
		const as3double beta  = rk4b[iStageRK];

    // Perform an entire SSPRK3 sweep and update the residual.
    UpdateTime(config_container,
               geometry_container,
               iteration_container,
               solver_container,
               spatial_container,
               initial_container,
               localTime, dtTime,
               alpha, beta,
               MonitoringData);
	}
}



void CSSPRK3Temporal::UpdateTime
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CIteration             *iteration_container,
 CSolver                *solver_container,
 CSpatial               *spatial_container,
 CInitial               *initial_container,
 as3double               localTime,
 as3double               dtTime,
 as3double               alpha,
 as3double               beta,
 as3vector1d<as3double> &MonitoringData

)
 /*
	* Function that performs a single SSPRK3 stage sweep and update the residual.
	*/
{
  // Initialize max of the Mach number squared.
  as3double M2max = 0.0;

  // Extract current total solution.
  auto& sol = solver_container->GetDataSolution();
  // Extract current total residual.
  auto& res = solver_container->GetDataResidual();
  // Extract tentative solution.
  auto& tmp = DataSolutionTentative;

  // Total solution size in first dimension.
  unsigned short nSolSize = sol.size();

  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

    // Loop over all nodes and copy solution to physical nodes of working array.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3)
#endif
    // Loop over every variable and every node.
    for(unsigned short iVar=0; iVar<nSolSize; iVar++){
      for(unsigned long j=0; j<nyNode; j++){
        for(unsigned long i=0; i<nxNode; i++){
          // Ghost node index.
          unsigned long I = i + NxStencil;
          unsigned long J = j + NyStencil;

          // Equivalent 1D running index of ghost node.
          unsigned long IJ = J*nxTotal + I;
          // Equivalent 1D running index of physical node.
          unsigned long ij = j*nxNode  + i;

          // Copy data point.
          work_array[iVar][IJ] = sol[iVar][ij];
        }
      }
    }


    // Impose the boundary conditions on all boundary nodes.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned short iBoundary=0; iBoundary<nFace; iBoundary++){

      // Extract the relevant boundary container.
      auto* boundary_container = solver_container->GetBoundaryContainer(iBoundary);

      // Impose boundary condition.
      boundary_container->ImposeBoundaryCondition(config_container,
                                                  geometry_container,
                                                  work_array,
                                                  localTime);
    }

  } // End of OpenMP parallel region.


  // Preprocess the data, if needed.
  iteration_container->Preprocess(config_container,
                                  geometry_container,
                                  solver_container,
                                  spatial_container,
                                  work_array,
                                  localTime);

  // Compute the residual on all physical nodes.
  spatial_container->ComputeResidual(config_container,
                                     geometry_container,
                                     spatial_container,
                                     initial_container,
                                     work_array,
                                     res,
                                     localTime,
                                     MonitoringData);

  // Set the max of the Mach squared in this element.
  M2max = MonitoringData[0];

  // Abbreviations for updating the residual according to a SSPRK3 scheme.
  const as3double oma = 1.0 - alpha;
  const as3double bdt = beta*dtTime;

  // Distinguish between the first stage evaluation and the rest.
  if( fabs( oma ) < 1.0e-10 ){

    // Loop over all nodes and update the solution.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(2)
#endif
    // Loop over every variable and node.
    for(unsigned short iVar=0; iVar<nSolSize; iVar++){
      for(unsigned long l=0; l<nNode; l++){
        // Predictor step, update the tentative solution.
        tmp[iVar][l]  = sol[iVar][l];
        // Corrector step, update the true solution.
        sol[iVar][l] += bdt*res[iVar][l];
      }
    } // End of OpenMP parallel region.
  }
  else {
    // This is the second or third stage of the SSPRK3.

    // Loop over all nodes and update the solution.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(2)
#endif
    // Loop over every variable and node.
    for(unsigned short iVar=0; iVar<nSolSize; iVar++){
      for(unsigned long l=0; l<nNode; l++){
        // Predictor-orrector step, update the true solution.
        sol[iVar][l] =   oma*sol[iVar][l]
                     + alpha*tmp[iVar][l]
                     +   bdt*res[iVar][l];
      }
    } // End of OpenMP parallel region.
  }


  // Assign the actual max of the Mach number.
  MonitoringData[0] = sqrt(M2max);
}