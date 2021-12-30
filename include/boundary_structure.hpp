#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"



class CBoundary {

	public:
		// Constructor.
		CBoundary(CConfig       *config_container,
							CGeometry     *geometry_container,
							unsigned short iBoundary);

		// Destructor.
		virtual ~CBoundary(void);

    // Pure virtual function that applies the boundary condition.
    // Must be overriden by a derived class.
    virtual void ImposeBoundaryCondition(CConfig              *config_container,
                                         CGeometry            *geometry_container,
                                         as3data1d<as3double> &work_array,
                                         as3double             localTime) = 0;

	protected:
		// Current boundary ID.
		unsigned short iBoundaryID;

	private:

};


class CEEBoundary : public CBoundary {

	public:
		// Constructor.
		CEEBoundary(CConfig       *config_container,
								CGeometry     *geometry_container,
								unsigned short iBoundary);

		// Destructor.
		virtual ~CEEBoundary(void) override;

  protected:

  private:

};


class CEEInterfaceBoundary : public CEEBoundary {

	public:
		// Constructor.
		CEEInterfaceBoundary(CConfig       *config_container,
												 CGeometry     *geometry_container,
												 unsigned short iBoundary);

		// Destructor.
		~CEEInterfaceBoundary(void) final;

    // Function that applies an interface boundary condition.
    void ImposeBoundaryCondition(CConfig              *config_container,
                                 CGeometry            *geometry_container,
                                 as3data1d<as3double> &work_array,
                                 as3double             localTime) final;
	protected:
		// Matching boundary ID.
		unsigned short jBoundaryID;

    // Nodal indices that pair the matching boundary.
    as3vector2d<unsigned long> NodeIndexPair;

	private:
    // Function that pairs the periodic data to the ghost nodes.
    void PairPeriodicIndices(CConfig   *config_container,
                             CGeometry *geometry_container);

    // Function that reports information on this boundary.
    void ReportOutput(void);
};


