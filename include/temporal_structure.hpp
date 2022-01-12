#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "iteration_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"



class CTemporal {

	public:
		// Constructor.
		CTemporal(CConfig     *config_container,
							CGeometry   *geometry_container,
							CIteration  *iteration_container,
							CSolver     *solver_container,
							CSpatial    *spatial_container);

		// Destructor.
		virtual ~CTemporal(void);

		// Pure virtual Function that performs an update over the entire simulation in time.
		// Must be overriden by one of the derived classes.
		virtual void TimeMarch(CConfig                *config_container,
        					 				 CGeometry              *geometry_container,
        					 				 CIteration             *iteration_container,
        					 				 CSolver                *solver_container,
        					 				 CSpatial               *spatial_container,
                           CInitial               *initial_container,
        					 				 as3double               physicalTime,
        									 as3double               dtTime,
                           as3vector1d<as3double> &MonitoringData) = 0;

	protected:
    // Number of nodes in x-direction.
    unsigned long nxNode;
    // Number of nodes in y-direction.
    unsigned long nyNode;
		// Total number of nodes.
		unsigned long nNode;

    // Total number of nodes in buffer layer.
    unsigned long nNodeBuffer;

    // Total nodes in x-direction, including ghost nodes.
    unsigned long nxTotal;
    // Total nodes in y-direction, including ghost nodes.
    unsigned long nyTotal;

    // Stencil size in west(left) x-direction.
		unsigned short NxStencil;
		// Stencil size in east(right) x-direction.
		unsigned short MxStencil;

		// Stencil size in south(bottom) y-direction.
		unsigned short NyStencil;
		// Stencil size in north(top) y-direction.
		unsigned short MyStencil;

    // Number of degrees-of-freedom in the working array.
    unsigned long  nWorkingArrayDOFs;
    // Number of variables needed in the working array.
    unsigned short nWorkingArrayVar;
    // Number of data entries needed in the working array.
    unsigned short nWorkingArrayEntry;

    // Work array needed for performing each grid sweep.
    as3data1d<as3double> work_array;

    // Function that initializes and defined the dimension parameters needed
    // in the working array.
    void InitializeWorkArrayDimension(CConfig *config_container);

    // Function that initializes the working array by reserving its memory.
    void InitializeWorkArray(void);

	private:

};



class CLSRK4Temporal : public CTemporal {

	public:
		// Constructor.
		CLSRK4Temporal(CConfig     *config_container,
    							 CGeometry   *geometry_container,
    							 CIteration  *iteration_container,
    							 CSolver     *solver_container,
    							 CSpatial    *spatial_container);

		// Destructor.
		~CLSRK4Temporal(void) final;

		// Function that performs an update over the entire simulation in time.
		void TimeMarch(CConfig                *config_container,
					 				 CGeometry              *geometry_container,
					 				 CIteration             *iteration_container,
					 				 CSolver                *solver_container,
					 				 CSpatial               *spatial_container,
                   CInitial               *initial_container,
					 				 as3double               physicalTime,
									 as3double               dtTime,
                   as3vector1d<as3double> &MonitoringData) final;

	protected:

	private:
		// Number of RK stages: this is a 5-stage scheme.
		unsigned short nStageRK = LSRK4_N_STAGES;
		// Number of storage for tentative data.
		unsigned short nStorage = LSRK4_N_STORAGE;

		// LSRK4 coefficients.
		as3vector1d<as3double> rk4a;
		as3vector1d<as3double> rk4b;
		as3vector1d<as3double> rk4c;

		// Tentative/intermediate physical solution stored.
		// Dimension: [iVar][iNode].
		as3data1d<as3double> DataSolutionTentative;

    // Tentative/intermediate auxiliary solution stored.
		// Dimension: [iVar][iNode].
		as3data1d<as3double> DataSolutionTentativeAux;

    // Function that performs a single stage sweep of a LSRK4.
		void UpdateTime(CConfig                *config_container,
									  CGeometry              *geometry_container,
									  CIteration             *iteration_container,
									  CSolver                *solver_container,
									  CSpatial               *spatial_container,
                    CInitial               *initial_container,
									  as3double               localTime,
                    as3double               dtTime,
                    as3double               alpha,
                    as3double               beta,
                    as3vector1d<as3double> &MonitoringData);
};
