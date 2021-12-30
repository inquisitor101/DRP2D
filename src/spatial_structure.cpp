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
  // Extract number of nodes in x-direction.
  nxNode = config_container->GetnxNode();
  // Extract number of nodes in y-direction.
  nyNode = config_container->GetnyNode();
  // Total number of nodes in grid.
  nNode  = geometry_container->GetnNode();

  // Extract stencil widths, per dimension.
  MxStencil = config_container->GetmxStencil();
  NxStencil = config_container->GetnxStencil();
  MyStencil = config_container->GetmyStencil();
  NyStencil = config_container->GetnyStencil();

  // Total nodes in x-diretcion.
  nxTotal = nxNode + MxStencil + NxStencil;
  // Total nodes in y-diretcion.
  nyTotal = nyNode + MyStencil + NyStencil;

  // Initialize and define stencil information.
  InitializeStencilStrategy(config_container,  geometry_container);

  // Identify required nodal indices.
  IdentifyRequiredNodalIndices(config_container, geometry_container);
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


void CSpatial::IdentifyRequiredNodalIndices
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that identifies the needed indices, given the stencils used.
	*/
{
  // Required indices.
  unsigned long nRequired = nxTotal*nyTotal
                          - NxStencil*NyStencil  // corner: 0
                          - MxStencil*NyStencil  // corner: 1
                          - NxStencil*MyStencil  // corner: 2
                          - MxStencil*MyStencil; // corner: 3

  // Reserve memory for the required nodal indices.
  RequiredIndices.resize(nRequired);

  // Corner: 0 indicial locations.
  unsigned long C0x = NxStencil,           C0y = NyStencil;
  unsigned long C1x = nxTotal-MxStencil-1, C1y = NyStencil;
  unsigned long C2x = NxStencil,           C2y = nyTotal-MyStencil-1;
  unsigned long C3x = nxTotal-MxStencil-1, C3y = nyTotal-MyStencil-1;

  // Loop over all nodes and identify the relevant ones.
  unsigned long k = 0;
  for(unsigned long j=0; j<nyTotal; j++){
    for(unsigned long i=0; i<nxTotal; i++){

      // Check if this is corner: 0.
      if( i<C0x && j<C0y ) continue;
      // Check if this is corner: 1.
      if( i>C1x && j<C1y ) continue;
      // Check if this is corner: 2.
      if( i<C2x && j>C2y ) continue;
      // Check if this is corner: 3.
      if( i>C3x && j>C3y ) continue;

      // Deduce the 1D running index.
      unsigned long l = j*nxTotal + i;

      // Book-keep global running 1D index.
      RequiredIndices[k++] = l;
    }
  }
}


void CSpatial::InitializeStencilStrategy
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the finite-difference stencils.
	*/
{
  // Reserve memory for each dimension of the stencil indices.
  ijStencil.resize(nDim);
  // Reserve memory for each dimension for the stencil coefficients.
  cStencil.resize(nDim);

  // Assemble stencil information.
  for(unsigned short iDim=0; iDim<nDim; iDim++){

    // Relevant grid resolution.
    const as3double ds     = geometry_container->GetGridResolution(iDim);
    // Stencil size in backward x-direction.
    const unsigned short n = config_container->GetnnStencil(iDim);
    // Stencil size in forward  x-direction.
    const unsigned short m = config_container->GetmmStencil(iDim);

    // Check what type of stencil we are dealing with.
    switch( config_container->GetTypeStencil(iDim) ){

      case(STENCIL_DRP_M3N3):
      {
        // Non-zero coefficients in stencil.
        unsigned short ns = 6;

        // Allocate memory for indices.
        ijStencil[iDim].resize(ns);
        // Allocate memory for coefficients.
        cStencil[iDim].resize(ns);

        // Explicitly populate the stencil coefficients.
        ijStencil[iDim][0] = -3; cStencil[iDim][0] = -0.020843142770;
        ijStencil[iDim][1] = -2; cStencil[iDim][1] =  0.166705904414;
        ijStencil[iDim][2] = -1; cStencil[iDim][2] = -0.770882380518;
        ijStencil[iDim][3] =  1; cStencil[iDim][3] =  0.770882380518;
        ijStencil[iDim][4] =  2; cStencil[iDim][4] = -0.166705904414;
        ijStencil[iDim][5] =  3; cStencil[iDim][5] =  0.020843142770;

        // Consistency check.
        assert( ns == (m + n) );

        break;
      }

      case(STENCIL_DRP_M2N4):
      {
        // Non-zero coefficients in stencil.
        unsigned short ns = 7;

        // Allocate memory for indices.
        ijStencil[iDim].resize(ns);
        // Allocate memory for coefficients.
        cStencil[iDim].resize(ns);

        // Explicitly populate the stencil coefficients.
        ijStencil[iDim][0] = -4; cStencil[iDim][0] =  0.0161404967150957;
        ijStencil[iDim][1] = -3; cStencil[iDim][1] = -0.1228212790198640;
        ijStencil[iDim][2] = -2; cStencil[iDim][2] =  0.4553322777062210;
        ijStencil[iDim][3] = -1; cStencil[iDim][3] = -1.2492595882614900;
        ijStencil[iDim][4] =  0; cStencil[iDim][4] =  0.5018904380193460;
        ijStencil[iDim][5] =  1; cStencil[iDim][5] =  0.4399321927296360;
        ijStencil[iDim][6] =  2; cStencil[iDim][6] = -0.0412145378889463;

        // Consistency check.
        assert( ns == (m + n + 1) );

        break;
      }

      default:
        Terminate("CSpatial::InitializeStencilStrategy", __FILE__, __LINE__,
                  "Type of stencil unknown,");
    }

    // Factor in the grid resolution in the stencil coefficients.
    for(unsigned short i=0; i<cStencil[iDim].size(); i++) cStencil[iDim][i] /= ds;
  }
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


void CEESpatial::ComputeResidual
(
 CConfig                *config_container,
 CGeometry              *geometry_container,
 CSpatial               *spatial_container,
 CInitial               *initial_container,
 as3data1d<as3double>   &work_array,
 as3data1d<as3double>   &residual,
 as3double               localTime,
 as3vector1d<as3double> &MonitoringData
)
 /*
	* Function that sweeps through the grid and updates the residual in an EE solver.
	*/
{
  // Abbreviation.
  const as3double gm1 = GAMMA - 1.0;
  // Initialize max Mach number.
  as3double M2max     = 0.0;

  // Abbreviation for stencils.
  const unsigned short nxs = cStencil[0].size();
  const unsigned short nys = cStencil[1].size();
  const long          *is  = ijStencil[0].data();
  const long          *js  = ijStencil[1].data();
  const as3double     *cx  = cStencil[0].data();
  const as3double     *cy  = cStencil[1].data();

  // Indices corresponding to physical region.
  unsigned long I0 = NxStencil, I1 = I0 + nxNode;
  unsigned long J0 = NyStencil, J1 = J0 + nyNode;

  // Assign the pointers for the working solution of the fluxes.
  // Note, the working array is in the same address as the x-flux.
  as3double **dVarDx = work_array.data();
  as3double **dVarDy = dVarDx + nVar;


  // Initiate OpenMP parallel region, if specified.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

    // Pre-compute the x- and y-fluxes on all nodes.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), reduction(max:M2max)
#endif
    // Loop over all nodes and compute the y-flux.
    for(unsigned long k=0; k<RequiredIndices.size(); k++){

      // Deduce the 1D running index.
      unsigned long l = RequiredIndices[k];

      // Copy conservative data.
      const as3double rho  = work_array[0][l];
      const as3double rhou = work_array[1][l];
      const as3double rhov = work_array[2][l];
      const as3double rhoE = work_array[3][l];

      // Compute primitive variables.
      const as3double ovrho = 1.0/rho;
      const as3double u     = ovrho*rhou;
      const as3double v     = ovrho*rhov;
      const as3double p     = gm1*( rhoE - 0.5*(u*rhou + v*rhov) );

      // Magnitude of the velocity squared.
      const as3double umag2 = u*u + v*v;
      // Speed of sound squared.
      const as3double a2    = GAMMA*p*ovrho;

      // Compute the local Mach number squared.
      const as3double M2 = umag2/a2;
      // Check if this value is the largest.
      M2max = std::max(M2max, M2);

      // Compute and store the y-flux.
      dVarDy[0][l] =     rhov;
      dVarDy[1][l] = u*  rhov;
      dVarDy[2][l] = v*  rhov + p;
      dVarDy[3][l] = v*( rhoE + p );

      // Compute and store the x-flux.
      dVarDx[0][l] =     rhou;
      dVarDx[1][l] = u*  rhou + p;
      dVarDx[2][l] = v*  rhou;
      dVarDx[3][l] = u*( rhoE + p );
    }


    // Compute the residual.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3)
#endif
    // Loop over all nodes and variables and update the residual.
    for(unsigned short iVar=0; iVar<nVar; iVar++){
      for(unsigned long j=J0; j<J1; j++){
        for(unsigned long i=I0; i<I1; i++){

          // Equivalent 1D running node index.
          const unsigned long ij = j*nxTotal + i;

          // Compute x-derivative approximation.
          as3double tmpx = 0.0;
          for(auto k=0; k<nxs; k++) tmpx -= cx[k]*dVarDx[iVar][ ij + is[k] ];

          // Compute y-derivative approximation.
          as3double tmpy = 0.0;
          for(auto k=0; k<nys; k++) tmpy -= cy[k]*dVarDy[iVar][ ij +js[k]*nxTotal ];

          // Deduce 1D running index.
          unsigned long l = (j-NyStencil)*nxNode + (i-NxStencil);

          // Add the contrinution to the residual.
          residual[iVar][l] = tmpx + tmpy;

          // Deduce the 1D running index.
          l++;
        }
      }
    }

  } // End of OpenMP parallel loop.

  // Assign the monitoring data.
  MonitoringData[0] = M2max;
}



