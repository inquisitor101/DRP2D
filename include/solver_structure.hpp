#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "spatial_structure.hpp"
#include "initial_structure.hpp"



class CSolver {

	public:
		// Constructor.
		CSolver(CConfig   *config_container,
						CGeometry *geometry_container,
            CInitial  *initial_container,
						CSpatial  *spatial_container);

		// Destructor.
		virtual ~CSolver(void);

    // Pure virtual function that sets the initial condition.
		// Must be overriden by one of the derived classes.
		virtual void InitializeSolution(CGeometry *geometry_container,
                                    CInitial  *initial_container,
                                    as3double  time) = 0;

    // Pure virtual function that computes the maximum stable time step per element.
    // Must be overriden by a derived class.
    virtual as3double ComputeTimeStep(void) = 0;

    // Getter: returns DataSolution.
    as3data1d<as3double> &GetDataSolution(void)          {return DataSolution;}

	protected:
    // Total number of nodes.
    unsigned long nNode;
    // Number of boundary conditions in this zone.
    unsigned short nBoundary;

    // Data of residual at points in grid.
    // Dimension: [iVar][iNode].
    as3data1d<as3double> DataResidual;

    // Data of solution points in grid.
    // Dimension: [iVar][iNode].
    as3data1d<as3double> DataSolution;

    // Data of fluxes as a function of the solution points in grid.
    // Dimension: [iDim][iVar][iNode].
    as3data2d<as3double> DataFlux;

	private:

};


class CEESolver : public CSolver {

	public:
		// Constructor.
		CEESolver(CConfig   *config_container,
							CGeometry *geometry_container,
              CInitial  *initial_container,
							CSpatial  *spatial_container);

		// Destructor.
		~CEESolver(void) final;

    // Function that sets the initial condition.
    void InitializeSolution(CGeometry *geometry_container,
                            CInitial  *initial_container,
                            as3double  time) final;

    // Function that computes the maximum stable time step per element.
    // Must be overriden by a derived class.
    as3double ComputeTimeStep(void);

	protected:

	private:

};