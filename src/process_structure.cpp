#include "process_structure.hpp"



CProcess::CProcess
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 COutput       *output_container,
 CInitial      *initial_container,
 CStencil     **stencil_container,
 CSolver       *solver_container,
 CSpatial      *spatial_container,
 unsigned short iZone
)
 /*
  * Constructor, used to initialize CProcess.
  */
{

}


CProcess::~CProcess
(
  void
)
 /*
  * Destructor for CProcess class, frees allocated memory.
  */
{

}


void CProcess::FilterSolution
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSolver   *solver_container,
 CSpatial  *spatial_container
)
 /*
  * Function that filters the solution.
  */
{

}


CEEProcess::CEEProcess
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 COutput       *output_container,
 CInitial      *initial_container,
 CStencil     **stencil_container,
 CSolver       *solver_container,
 CSpatial      *spatial_container,
 unsigned short iZone
)
  :
    CProcess
    (
      config_container,
      geometry_container,
      output_container,
      initial_container,
      stencil_container,
      solver_container,
      spatial_container,
      iZone
    )
 /*
  * Constructor, used to initialize CEEProcess.
  */
{

}


CEEProcess::~CEEProcess
(
  void
)
 /*
  * Destructor for CEEProcess class, frees allocated memory.
  */
{
  for(unsigned short i=0; i<sensor_container.size(); i++)
    if( sensor_container[i] ) delete sensor_container[i];
}


void CEEProcess::WriteProcessedData
(
  CConfig   *config_container,
  CGeometry *geometry_container,
  COutput   *output_container
)
 /*
  * Function that writes the out the processed data to a file.
  */
{

}


void CEEProcess::ProcessData
(
  CConfig   *config_container,
  CGeometry *geometry_container,
  CSpatial  *spatial_container,
  CSolver   *solver_container,
  CInitial  *initial_container,
  as3double  localTime
)
 /*
  * Function that selects what/how to process the data.
  */
{

}


void CEEProcess::ComputeReflectionCoefficientFuncP
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that computes the reflection coefficient as a function of (P).
  * That is: R = max( abs( p'(t)-p(t=0) ) )/p'(t=0).
  */
{

}


void CEEProcess::ComputeReflectionCoefficientFuncPT
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that computes the reflection coefficient as a function of (P, T).
  * That is: R = ( max( p'(t) )/pinf ) / ( max( T'(t=0) )/ Tinf ).
  */
{

}


void CEEProcess::ComputeReflectionCoefficientFuncPM
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that computes the reflection coefficient as a function of (P, M).
  * That is: R = ( max( p(t)/pinf ) ) / ( max( M(t) )/ Minf ).
  */
{

}


void CEEProcess::SampleProbePoints
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSpatial  *spatial_container,
 CSolver   *solver_container,
 CInitial  *initial_container,
 as3double  localTime
)
 /*
  * Function that samples the pressure and v-velocity at different probe locations.
  */
{

}


CSensorData::CSensorData
(
 CConfig               *config_container,
 CGeometry             *geometry_container,
 as3vector1d<as3double> probe
)
 /*
  * Constructor used to initialize CSensorData.
  */
{

}


CSensorData::~CSensorData
(
 void
)
 /*
  * Destructor for CSensorData class, frees allocated memory.
  */
{

}



