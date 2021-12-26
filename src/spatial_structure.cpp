#include "spatial_structure.hpp"



CSpatial::CSpatial
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CInitial  *initial_container
)
 /*
	* Constructor, used to initialize CSpatial.
	*/
{

}


CSpatial::~CSpatial
(
 void
)
 /*
	* Destructor for CSpatial class, frees allocated memory.
	*/
{

}


CEESpatial::CEESpatial
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CInitial  *initial_container
)
	:
		CSpatial
		(
		 config_container,
		 geometry_container,
     initial_container
		)
 /*
	* Constructor, used to initialize CEESpatial.
	*/
{

}


CEESpatial::~CEESpatial
(
 void
)
 /*
	* Destructor for CEESpatial class, frees allocated memory.
	*/
{

}
