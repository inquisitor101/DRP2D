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
  // Total number of nodes in grid.
  nNode = geometry_container->GetnNode();
}


CTemporal::~CTemporal
(
 void
)
 /*
	* Destructor for CTemporal class, frees allocated memory.
	*/
{

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
}


CLSRK4Temporal::~CLSRK4Temporal
(
 void
)
 /*
	* Destructor for CLSRK4Temporal class, frees allocated memory.
	*/
{

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
//   // Initialize max of the Mach number squared.
//   as3double M2max = 0.0;
//
//   // Initiate OpenMP parallel region, if specified.
// #ifdef HAVE_OPENMP
// #pragma omp parallel
// #endif
//   {
//     // Work array needed for performing a single grid sweep.
//     as3data1d<as3double> work_array;
//     // Reserve needed memory for the working array.
//     InitializeWorkArray(work_array);
//
//
//   // Impose the boundary conditions in all boundary elements across all zones.
// #ifdef HAVE_OPENMP
// #pragma omp for schedule(static), collapse(2)
// #endif
//     for(unsigned short iZone=0; iZone<nZone; iZone++){
//       for(unsigned short iBoundary=0; iBoundary<nFace; iBoundary++){
//
//         // Extract the relevant boundary container.
//         auto* boundary_container = solver_container[iZone]->GetBoundaryContainer(iBoundary);
//
//         // Impose boundary condition.
//         boundary_container->ImposeBoundaryCondition(config_container,
//                                                     geometry_container,
//                                                     solver_container,
//                                                     element_container,
//                                                     spatial_container,
//                                                     localTime);
//       }
//     }
//
//
//     // Loop over all the elements in all the zones and compute their residual.
// #ifdef HAVE_OPENMP
// #pragma omp for schedule(static), reduction(max:M2max)
// #endif
//     for(unsigned long i=0; i<nElemTotal; i++){
//
//       // Extract local zone number.
//       unsigned short iZone = MapGlobalToLocal[i][0];
//       // Extract local zone element index.
//       unsigned long iElem  = MapGlobalToLocal[i][1];
//
//       // Preprocess the data and apply the boundary conditions.
//       iteration_container[iZone]->Preprocess(config_container,
//   																					 geometry_container,
//   																					 solver_container,
//   																					 element_container,
//   																					 spatial_container,
//                                              work_array,
//   																					 localTime, iElem);
//
//       // Perform a single grid sweep over a single unique element.
//       iteration_container[iZone]->ComputeResidual(config_container,
//       																						geometry_container,
//       																						solver_container[iZone],
//       																						element_container[iZone],
//       																						spatial_container[iZone],
//                                                   initial_container[iZone],
//                                                   work_array,
//       																						localTime, iElem,
//                                                   MonitoringData);
//
//       // Set the max of the Mach squared in this element.
//       M2max = std::max(M2max, MonitoringData[0]);
//     }
//
//
//     // Loop over all elements in all zones and update the solution.
// #ifdef HAVE_OPENMP
// #pragma omp for schedule(static)
// #endif
//     for(unsigned long i=0; i<nElemTotal; i++){
//
//       // Extract local zone number.
//       unsigned short iZone = MapGlobalToLocal[i][0];
//       // Extract local zone element index.
//       unsigned long iElem  = MapGlobalToLocal[i][1];
//
//       // Extract current data container.
//       auto* data_container = solver_container[iZone]->GetDataContainer(iElem);
//
//       // Number of nodes per zone.
//       unsigned short nNode = nNodeZone[iZone];
//
//       // Extract current total solution.
//       auto& sol = data_container->GetDataDOFsSol();
//       // Extract current total residual.
//       auto& res = data_container->GetDataDOFsRes();
//       // Extract tentative solution.
//       auto& tmp = DataDOFsSolTentative[iZone][iElem];
//
//       // Loop over every variable.
//       for(unsigned short iVar=0; iVar<sol.size(); iVar++){
//
//         // Loop over nodes and update the residual.
// #pragma omp simd
//         for(unsigned short iNode=0; iNode<nNode; iNode++){
//           // Predictor step, update the tentative solution.
//           tmp[iVar][iNode]  = alpha*tmp[iVar][iNode] + dtTime*res[iVar][iNode];
//         }
// #pragma omp simd
//         for(unsigned short iNode=0; iNode<nNode; iNode++){
//           // Corrector step, update the true solution.
//           sol[iVar][iNode] += beta*tmp[iVar][iNode];
//         }
//       }
//     }
//
//     // Free the local work array.
//     for(unsigned short i=0; i<work_array.size(); i++)
//       if( work_array[i] ) delete [] work_array[i];
//
//   } // End of OpenMP parallel region.
//
//   // Assign the actual max of the Mach number.
//   MonitoringData[0] = sqrt(M2max);
}

