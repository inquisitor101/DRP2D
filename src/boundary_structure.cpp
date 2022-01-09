#include "boundary_structure.hpp"



CBoundary::CBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CStencil     **stencil_container,
 unsigned short iBoundary
)
 /*
	* Constructor, used to initialize CBoundary.
	*/
{
	// Assign boundary ID that is unique to this grid.
	iBoundaryID = iBoundary;
}


CBoundary::~CBoundary
(
 void
)
 /*
	* Destructor for CBoundary class, frees allocated memory.
	*/
{

}


CEEBoundary::CEEBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CStencil     **stencil_container,
 unsigned short iBoundary
)
	:
		CBoundary
		(
		 config_container,
		 geometry_container,
     stencil_container,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEBoundary.
	*/
{

}


CEEBoundary::~CEEBoundary
(
 void
)
 /*
	* Destructor for CEEBoundary class, frees allocated memory.
	*/
{

}


CEEInterfaceBoundary::CEEInterfaceBoundary
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 CStencil     **stencil_container,
 unsigned short iBoundary
)
	:
		CEEBoundary
		(
		 config_container,
		 geometry_container,
     stencil_container,
		 iBoundary
		)
 /*
	* Constructor, used to initialize CEEInterfaceBoundary.
	*/
{
  // // Determine matching boundary.
  // switch(iBoundary){
  //   case IDX_SOUTH: jBoundaryID = IDX_NORTH; break;
  //   case IDX_NORTH: jBoundaryID = IDX_SOUTH; break;
  //   case IDX_WEST:  jBoundaryID = IDX_EAST;  break;
  //   case IDX_EAST:  jBoundaryID = IDX_WEST;  break;
  //   default:
  //     Terminate("CEEInterfaceBoundary::CEEInterfaceBoundary", __FILE__, __LINE__,
  //               "Unknown matching boundary ID.");
  // }
  //
  // // Pair and map between the periodic and ghost indices.
  // PairPeriodicIndices(config_container,
  //                     geometry_container);
  //
  // // Report output.
  // ReportOutput();
}


CEEInterfaceBoundary::~CEEInterfaceBoundary
(
 void
)
 /*
	* Destructor for CEEInterfaceBoundary class, frees allocated memory.
	*/
{

}


void CEEInterfaceBoundary::ImposeBoundaryCondition
(
 CConfig              *config_container,
 CGeometry            *geometry_container,
 as3data1d<as3double> &work_array,
 as3double             localTime
)
 /*
	* Function that imposes the current interface boundary condition.
	*/
{
//   // Loop over all the relevant nodal data and copy the paired periodic boundary stencil.
//   for(unsigned short iVar=0; iVar<work_array.size(); iVar++){
// #pragma omp simd
//     for(unsigned long l=0; l<NodeIndexPair.size(); l++){
//       work_array[iVar][ NodeIndexPair[l][0] ] = work_array[iVar][ NodeIndexPair[l][1] ];
//     }
//   }
}


void CEEInterfaceBoundary::ReportOutput
(
 void
)
 /*
  * Function that reports output for monitoring.
  */
{
  // Report output.
  std::cout << "  "
            << "iBoundary(" << iBoundaryID      << "): "
            << DisplayBoundarySide(iBoundaryID) << " "
            << "jBoundary(" << jBoundaryID      << "): "
            << DisplayBoundarySide(jBoundaryID)
            << std::endl;
}


void CEEInterfaceBoundary::PairPeriodicIndices
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
  * Function that pairs the periodic data to the ghost nodes.
  */
{
  // // Extract number of nodes in x-direction.
  // unsigned long nxNode = config_container->GetnxNode();
  // // Extract number of nodes in y-direction.
  // unsigned long nyNode = config_container->GetnyNode();
  //
  // // Extract stencil widths, per dimension.
  // unsigned short Mx = config_container->GetmxStencil();
  // unsigned short Nx = config_container->GetnxStencil();
  // unsigned short My = config_container->GetmyStencil();
  // unsigned short Ny = config_container->GetnyStencil();
  //
  // // Total nodal points per dimension, including ghost nodes.
  // unsigned long nxTotal = nxNode + Mx + Nx;
  // unsigned long nyTotal = nyNode + My + Ny;
  //
  // // Temporary storage of dimensionally-split nodal indices, for readability.
  // as3vector2d<unsigned long> tmp;
  //
  // // Check which boundary indices and stencil to use.
  // switch(iBoundaryID){
  //
  //   // This is a south boundary with ghost nodes that need to be linked to the
  //   // north boundary stencil.
  //   case(IDX_SOUTH):
  //   {
  //     // Allocate memory for south-most nodal (ghost) indices.
  //     tmp.resize(nxNode*Ny, as3vector1d<unsigned long> (4) );
  //
  //     // Indices of the ghost nodes.
  //     unsigned long I0 = Nx;
  //     unsigned long I1 = nxTotal - Mx;
  //
  //     // Pair ghost-to-physical nodes.
  //     unsigned long k = 0;
  //     for(unsigned long jj=0; jj<Ny; jj++){
  //       for(unsigned long ii=I0; ii<I1; ii++){
  //
  //         // Indices of ghost nodes.
  //         tmp[k][0] = ii;
  //         tmp[k][1] = Ny - 1 - jj;
  //         // Indices of physical nodes.
  //         tmp[k][2] = ii;
  //         tmp[k][3] = nyTotal - 1 - My - jj;
  //
  //         // Increment global index.
  //         k++;
  //       }
  //     }
  //
  //     break;
  //   }
  //
  //   // This is a north boundary with ghost nodes that need to be linked to the
  //   // south boundary stencil.
  //   case(IDX_NORTH):
  //   {
  //     // Allocate memory for north-most nodal (ghost) indices.
  //     tmp.resize(nxNode*My, as3vector1d<unsigned long> (4) );
  //
  //     // Indices of the ghost nodes.
  //     unsigned long I0 = Nx;
  //     unsigned long I1 = nxTotal - Mx;
  //
  //     // Pair ghost-to-physical nodes.
  //     unsigned long k = 0;
  //     for(unsigned long jj=0; jj<My; jj++){
  //       for(unsigned long ii=I0; ii<I1; ii++){
  //
  //         // Indices of ghost nodes.
  //         tmp[k][0] = ii;
  //         tmp[k][1] = nyTotal - My + jj;
  //         // Indices of physical nodes.
  //         tmp[k][2] = ii;
  //         tmp[k][3] = Ny + jj;
  //
  //         // Increment global index.
  //         k++;
  //       }
  //     }
  //
  //     break;
  //   }
  //
  //   // This is a west boundary with ghost nodes that need to be linked to the
  //   // east boundary stencil.
  //   case(IDX_WEST):
  //   {
  //     // Allocate memory for west-most nodal (ghost) indices.
  //     tmp.resize(nyNode*Nx, as3vector1d<unsigned long> (4) );
  //
  //     // Indices of the ghost nodes.
  //     unsigned long J0 = Ny;
  //     unsigned long J1 = nyTotal - My;
  //
  //     // Pair ghost-to-physical nodes.
  //     unsigned long k = 0;
  //     for(unsigned long jj=J0; jj<J1; jj++){
  //       for(unsigned long ii=0; ii<Nx; ii++){
  //
  //         // Indices of ghost nodes.
  //         tmp[k][0] = Nx - 1 - ii;
  //         tmp[k][1] = jj;
  //         // Indices of physical nodes.
  //         tmp[k][2] = nxTotal - 1 - Mx - ii;
  //         tmp[k][3] = jj;
  //
  //         // Increment global index.
  //         k++;
  //       }
  //     }
  //
  //     break;
  //   }
  //
  //   // This is an east boundary with ghost nodes that need to be linked to the
  //   // west boundary stencil.
  //   case(IDX_EAST):
  //   {
  //     // Allocate memory for east-most nodal (ghost) indices.
  //     tmp.resize(nyNode*Mx, as3vector1d<unsigned long> (4) );
  //
  //     // Indices of the ghost nodes.
  //     unsigned long J0 = Ny;
  //     unsigned long J1 = nyTotal - My;
  //
  //     // Pair ghost-to-physical nodes.
  //     unsigned long k = 0;
  //     for(unsigned long jj=J0; jj<J1; jj++){
  //       for(unsigned long ii=0; ii<Mx; ii++){
  //
  //         // Indices of ghost nodes.
  //         tmp[k][0] = nxTotal - Mx + ii;
  //         tmp[k][1] = jj;
  //         // Indices of physical nodes.
  //         tmp[k][2] = Nx + ii;
  //         tmp[k][3] = jj;
  //
  //         // Increment global index.
  //         k++;
  //       }
  //     }
  //
  //     break;
  //   }
  //
  //   default:
  //     Terminate("CEEInterfaceBoundary::PairPeriodicIndices", __FILE__, __LINE__,
  //               "Type of boundary and stencil cannot be periodically paired.");
  // }
  //
  // // Convert the dimensionally-split nodal indices to a 1D running index.
  // NodeIndexPair.resize( tmp.size(), as3vector1d<unsigned long> (2) );
  //
  // // Populate the running indices of the ghost nodes.
  // for(unsigned long l=0; l<NodeIndexPair.size(); l++){
  //
  //   // Extract i-index of ghost node.
  //   unsigned long I = tmp[l][0];
  //   // Extract j-index of ghost node.
  //   unsigned long J = tmp[l][1];
  //   // Convert dimensionally-split ghost index to a global running index.
  //   NodeIndexPair[l][0] = J*nxTotal + I;
  //
  //   // Extract i-index of physical node.
  //   unsigned long i = tmp[l][2];
  //   // Extract j-index of physical node.
  //   unsigned long j = tmp[l][3];
  //   // Convert dimensionally-split physical index to a global running index.
  //   NodeIndexPair[l][1] = j*nxTotal + i;
  // }
}


