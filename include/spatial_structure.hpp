#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "initial_structure.hpp"


class CSpatial {

	public:
		// Constructor.
		CSpatial(CConfig   	*config_container,
						 CGeometry 	*geometry_container,
             CInitial   *initial_container);

		// Destructor.
		virtual ~CSpatial(void);

};


class CEESpatial : public CSpatial {

	public:
		// Constructor.
		CEESpatial(CConfig   	*config_container,
						 	 CGeometry 	*geometry_container,
               CInitial   *initial_container);

		// Destructor.
		~CEESpatial(void);

};



