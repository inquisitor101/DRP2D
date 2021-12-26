#pragma once

#include "option_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"



class CInput {

	public:
		// Constructor.
		CInput(CConfig   *config_container,
					 CGeometry *geometry_container);

		// Destructor.
		~CInput(void);

    // Function that reads the solution from a restart file.
    void ReadSolutionRestartFile(CConfig    *config_container,
                                 CGeometry  *geometry_container,
                                 CSolver    *solver_container,
                                 as3double  &SimTime);

	protected:

	private:
    
};



