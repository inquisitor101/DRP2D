#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"



class CIteration {

	public:
		// Constructor.
		CIteration(CConfig   *config_container,
							 CGeometry *geometry_container,
							 CSolver   *solver_container);

		// Destructor.
		virtual ~CIteration(void);

		// Pure virtual function that preprocesses every iteration.
		// Must be overriden by a derived class.
		virtual void Preprocess(CConfig              *config_container,
														CGeometry            *geometry_container,
														CSolver              *solver_container,
														CSpatial             *spatial_container,
                            as3data1d<as3double> &work_array,
														as3double             localTime) = 0;

    // Pure virtual function that sweeps through grid and updates residual accordingly.
    // Must be overriden by a derived class.
    virtual void ComputeResidual(CConfig                *config_container,
                                 CGeometry              *geometry_container,
                                 CSolver                *solver_container,
                                 CSpatial               *spatial_container,
                                 CInitial               *initial_container,
                                 as3data1d<as3double>   &work_array,
                                 as3double               localTime,
                                 as3vector1d<as3double> &MonitoringData) = 0;

	protected:

	private:

};


class CEEIteration : public CIteration {

	public:
		// Constructor.
		CEEIteration(CConfig   *config_container,
							 	 CGeometry *geometry_container,
							 	 CSolver   *solver_container);

		// Destructor.
		~CEEIteration(void) final;

	protected:
		// Function that preprocesses every iteration.
		void Preprocess(CConfig              *config_container,
										CGeometry            *geometry_container,
										CSolver              *solver_container,
										CSpatial             *spatial_container,
                    as3data1d<as3double> &work_array,
										as3double             localTime) final;

    // Function that sweeps through grid and updates residual accordingly.
    void ComputeResidual(CConfig                *config_container,
                         CGeometry              *geometry_container,
                         CSolver                *solver_container,
                         CSpatial               *spatial_container,
                         CInitial               *initial_container,
                         as3data1d<as3double>   &work_array,
                         as3double               localTime,
                         as3vector1d<as3double> &MonitoringData) final;

	private:

};


