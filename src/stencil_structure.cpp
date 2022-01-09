#include "stencil_structure.hpp"



CStencil::CStencil
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor, used to initialize CStencil class.
	*/
{

}


CStencil::~CStencil
(
 void
)
 /*
	* Deconstructor for CStencil class.
	*/
{

}


CDRPM3N3Stencil::CDRPM3N3Stencil
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
  :
    CStencil
    (
     config_container,
     geometry_container
    )
 /*
	* Constructor, used to initialize CDRPM3N3Stencil class.
	*/
{
  // Number of non-zero stencil in backward direction.
  NNStencil = 3;
  // Number of non-zero stencil in forward  direction.
  MMStencil = 3;

  // Stencil coefficients.
  cStencil.resize( MMStencil + NNStencil );
  // Stencil offset indices.
  kStencil.resize( MMStencil + NNStencil );

  // Populate stencil characteristics.
  kStencil[0] = -3; cStencil[0] = -0.020843142770;
  kStencil[1] = -2; cStencil[1] =  0.166705904414;
  kStencil[2] = -1; cStencil[2] = -0.770882380518;
  kStencil[3] =  1; cStencil[3] =  0.770882380518;
  kStencil[4] =  2; cStencil[4] = -0.166705904414;
  kStencil[5] =  3; cStencil[5] =  0.020843142770;
}


CDRPM3N3Stencil::~CDRPM3N3Stencil
(
 void
)
 /*
	* Deconstructor for CDRPM3N3Stencil class.
	*/
{

}











