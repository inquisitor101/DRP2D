#include "initial_structure.hpp"





CInitial::CInitial
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor, used to initialize CInitial.
	*/
{

}


CInitial::~CInitial
(
 void
)
 /*
	* Destructor for CInitial class, frees allocated memory.
	*/
{

}


CGaussianInitial::CGaussianInitial
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		CInitial
		(
		 config_container,
		 geometry_container
		)
 /*
	* Constructor, used to initialize CGaussianInitial.
	*/
{
  // Initialize Gaussian pulse center.
  x0 = config_container->GetCenterX0()[0];
  y0 = config_container->GetCenterX0()[1];
  // Initialize pulse strength, perctantage of background condition.
  A0 = config_container->GetDisturbanceRatio();
  // Initialize pulse radial distance.
  b  = config_container->GetDisturbanceWidth();
  // Pulse attenuation factor.
  kappa = log(2.0)/(b*b);

  // Background flow properties.
  Mach   = config_container->GetMachInf();
  Tinf   = 300.0;
  pInf   = 101325.0;
  // Deduce density.
  rhoInf = pInf/(GAS_CONSTANT*Tinf);

  // Reference speed of sound.
  aInf = sqrt(pInf*GAMMA/rhoInf);

  // Flow direction [degrees].
  theta  = config_container->GetFlowAngle();
  // Flow velocity.
  VelInf = Mach*aInf;
  uInf   = VelInf*cos(theta*PI_CONSTANT/180.0);
  vInf   = VelInf*sin(theta*PI_CONSTANT/180.0);

  // Initialize dispersion-relaxation correction.
  betaPML = uInf/( aInf*aInf - uInf*uInf );
}


CGaussianInitial::~CGaussianInitial
(
 void
)
 /*
	* Destructor for CGaussianInitial class, frees allocated memory.
	*/
{

}


as3vector1d<as3double> CGaussianInitial::SetInitialConditionPrimitive
(
 as3double x,
 as3double y,
 as3double t
)
 /*
	* Function that sets a Gaussian initial condition in this zone.
	*/
{
	// Compute relative position w.r.t. pulse center.
	as3double rxPos2 = (x-x0); rxPos2 *= rxPos2;
	as3double ryPos2 = (y-y0); ryPos2 *= ryPos2;

	// Radial distance squared.
	const as3double r2 = rxPos2 + ryPos2;

	// Compute initial condition in primitive form.
	const as3double rho = rhoInf;
	const as3double u   = uInf;
	const as3double v   = vInf;
	const as3double p   = pInf*( 1.0 + A0*exp(-kappa*r2) );

  // Return primitive-form vector.
  return {rho, u, v, p};
}


CIsentropicVortexInitial::CIsentropicVortexInitial
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		CInitial
		(
		 config_container,
		 geometry_container
		)
 /*
	* Constructor, used to initialize CIsentropicVortexInitial.
	*/
{
  // Initialize vortex center.
  x0 = config_container->GetCenterX0()[0];
  y0 = config_container->GetCenterX0()[1];
  // Initialize vortex.
  A0 = config_container->GetDisturbanceRatio();
  // Initialize vortex radial distance.
  Rv = config_container->GetDisturbanceWidth();

  // Background flow properties.
  Mach   = config_container->GetMachInf();
  Tinf   = 300.0;
  pInf   = 101325.0;
  // Deduce density.
  rhoInf = pInf/(GAS_CONSTANT*Tinf);

  // Reference speed of sound.
  aInf = sqrt(pInf*GAMMA/rhoInf);

  // Abbreviations.
  ovRv2 =  1.0/(Rv*Rv);
  alpha =  A0/(aInf*Rv);
  omega = -0.5*GAMMA*alpha*alpha;

  // Flow direction [degrees].
  theta  = config_container->GetFlowAngle();
  // Flow velocity.
  VelInf = Mach*aInf;
  uInf   = VelInf*cos(theta*PI_CONSTANT/180.0);
  vInf   = VelInf*sin(theta*PI_CONSTANT/180.0);

  // Initialize dispersion-relaxation correction.
  betaPML = uInf/( aInf*aInf - uInf*uInf );
}


CIsentropicVortexInitial::~CIsentropicVortexInitial
(
 void
)
 /*
	* Destructor for CIsentropicVortexInitial class, frees allocated memory.
	*/
{

}


as3vector1d<as3double> CIsentropicVortexInitial::SetInitialConditionPrimitive
(
 as3double x,
 as3double y,
 as3double t
)
 /*
	* Function that sets an isentropic vortex initial condition in this zone.
	*/
{
	// Compute relative position w.r.t. vortex center.
	const as3double dx  = x-x0; const as3double dx2 = dx*dx;
	const as3double dy  = y-y0; const as3double dy2 = dy*dy;

	// Radial distance squared.
	const as3double r2 = dx2 + dy2;

  // Compute stream function.
  const as3double psi = A0*exp(-0.5*r2*ovRv2);

	// Compute initial condition in primitive form.
	const as3double u   = -ovRv2*dy*psi + uInf;
	const as3double v   =  ovRv2*dx*psi + vInf;
	const as3double p   =  pInf*exp( omega*exp(-r2*ovRv2) );
  const as3double rho =  p/(GAS_CONSTANT*Tinf);

  // Return primitive-form vector.
  return {rho, u, v, p};
}

