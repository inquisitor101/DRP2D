#include "iteration_structure.hpp"



CIteration::CIteration
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSolver   *solver_container
)
 /*
	* Constructor, used to initialize CIteration.
	*/
{

}


CIteration::~CIteration
(
 void
)
 /*
	* Destructor for CIteration class, frees allocated memory.
	*/
{

}


CEEIteration::CEEIteration
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSolver   *solver_container
)
	:
		CIteration
		(
		 config_container,
		 geometry_container,
		 solver_container
		)
 /*
	* Constructor, used to initialize CEEIteration.
	*/
{

}


CEEIteration::~CEEIteration
(
 void
)
 /*
	* Destructor for CEEIteration class, frees allocated memory.
	*/
{

}


void CEEIteration::Preprocess
(
 CConfig              *config_container,
 CGeometry            *geometry_container,
 CSolver              *solver_container,
 CSpatial             *spatial_container,
 as3data1d<as3double> &work_array,
 as3double             localTime
)
 /*
	* Function that performs a preprocessing step for every iteration of the
  * EE solver.
	*/
{

}


void CEEIteration::ComputeResidual
(
 CConfig               *config_container,
 CGeometry             *geometry_container,
 CSolver               *solver_container,
 CSpatial              *spatial_container,
 CInitial              *initial_container,
 as3data1d<as3double>  &work_array,
 as3double              localTime,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that sweeps across space to update the residual in an EE solver.
	*/
{

}



