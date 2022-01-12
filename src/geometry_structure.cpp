#include "geometry_structure.hpp"



CGeometry::CGeometry
(
 CConfig *config_container
)
 /*
  * Constructor, Reads input grid.
  */
{
  // Extract number of nodes in x-direction.
  nxNode = config_container->GetnxNode();
  // Extract number of nodes in y-direction.
  nyNode = config_container->GetnyNode();
	// Deduce total number of nodes.
	nNode  = nxNode*nyNode;

  // Obtain grid domain size.
  GridSize = config_container->GetDomainBound();

  // Compute the total number of nPointCellP1.
  ComputePointCellP1();

  // Create the grid.
  GenerateGrid(nxNode, nyNode, nNode);

  // Generate buffer layer.
  GenerateBufferLayer(config_container);
}


CGeometry::~CGeometry
(
 void
)
 /*
  * Destructor for geometry class, frees allocated memory.
  */
{
  for(unsigned short i=0; i<GridCoordinate.size(); i++)
    if( GridCoordinate[i] ) delete [] GridCoordinate[i];
}


void CGeometry::ComputePointCellP1
(
 void
)
 /*
  * Function that computes the total number of points for a nPoly=1 cell-based grid.
  */
{
  // Total number of points needed to assemble an nPoly=1 cell grid.
  nPointCellP1 = (nxNode-1)*(nyNode-1)*N_POINTS_QUADRILATERAL;
}


void CGeometry::GenerateBufferLayer
(
 CConfig *config_container
)
 /*
  * Function that generates a buffer layer.
  */
{
  // Extract the starting location of the buffer layer in x-direction.
  const as3double xi = config_container->GetxBufferInterface();
  // Extract starting and ending location of the overall domain.
  const as3double x0 = GridSize[0], x1 = GridSize[1];

  // Consistency check.
  if( (xi >= x1) || (xi <= x0) )
    Terminate("CGeometry::GenerateBufferLayer", __FILE__, __LINE__,
              "Buffer layer position must be inside of the domain dimensions.");

  // Cast 1D array into 2D sub-arrays for readability, using the column-index.
  // That is: xcoord[nyNode][nxNode]
  // ... for efficiency (cpp row-major), loop over [nxNode] first.
  const as3double (*xcoord)[nxNode] = (const as3double (*)[nxNode]) GridCoordinate[0];

  // Loop over the grid and determine the i-direction buffer-layer indices.
  for(unsigned long i=0; i<nxNode; i++) if( xcoord[0][i] > xi ) IBufferLayer.emplace_back(i);

  // Number of nodes in buffer layer in x-direction.
  nxNodeBuffer = IBufferLayer.size();
  // Number of nodes in buffer layer in y-direction.
  nyNodeBuffer = nyNode;
  // Total number of nodes in the buffer layer.
  nNodeBuffer  = nxNodeBuffer*nyNodeBuffer;

  // Consistency check.
  if( !IBufferLayer.size() )
    Terminate("CGeometry::GenerateBufferLayer", __FILE__, __LINE__,
              "No nodes exist in buffer layer, something is wrong.");

  // Obtain damping function constant.
  const as3double sigma = config_container->GetDampingConstant();
  // Obtain damping function exponent.
  const as3double alpha = config_container->GetDampingExponent();

  // Reserve memory for the damping function along the x-direction.
  DampingFunctionXDir.resize(nxNodeBuffer, 0.0);

  // First half of the buffer layer.
  unsigned long ib = nxNodeBuffer/2;

  // Determine the inverse of the length of the first half of the buffer layer.
  const as3double ovDx1 = 1.0/( xcoord[0][ IBufferLayer[ib-1] ] - xcoord[0][ IBufferLayer[0]  ] );
  // Determine the inverse of the length of the second half of the buffer layer.
  const as3double ovDx2 = 1.0/( xcoord[0][ IBufferLayer.back()] - xcoord[0][ IBufferLayer[ib] ] );

  // Compute the damping function of the first half of the buffer layer.
  for(unsigned long i=0; i<ib; i++)
    DampingFunctionXDir[i] = sigma*pow( ovDx1*(xcoord[0][IBufferLayer[i]] - xcoord[0][IBufferLayer[0]]), alpha ) ;

  // Compute the damping function of the second half of the buffer layer.
  unsigned long k=1;
  for(unsigned long i=ib; i<nxNodeBuffer; i++)
    DampingFunctionXDir[nxNodeBuffer-k++] = sigma*pow( ovDx2*(xcoord[0][IBufferLayer[i]] - xcoord[0][IBufferLayer[ib]]), alpha ) ;
}


void CGeometry::GenerateGrid
(
 unsigned long nxNode,
 unsigned long nyNode,
 unsigned long nNode
)
 /*
  * Function that generates a grid.
  */
{
  // Reserve memory for coordinates.
  GridCoordinate.resize(nDim, nullptr);

  // Allocate in every variable.
  for(unsigned short iDim=0; iDim<GridCoordinate.size(); iDim++){

    // Allocate actual memory.
    GridCoordinate[iDim]  = new as3double[nNode]();

    // Check if allocation failed.
    if( !GridCoordinate[iDim] )
      Terminate("CGeometry::GenerateGrid", __FILE__, __LINE__,
                "Allocation failed for GridCoordinate.");
  }

  // Explicitly define grid boundary coordinates.
  const as3double x0 = GridSize[0], x1 = GridSize[1];
  const as3double y0 = GridSize[2], y1 = GridSize[3];

  // Compute the grid resolution.
  const as3double dx = (x1 - x0)/(nxNode-1.0);
  const as3double dy = (y1 - y0)/(nyNode-1.0);

  // Book-keep resolution.
  GridResolution[0] = dx;
  GridResolution[1] = dy;

  // Cast 1D array into 2D sub-arrays for readability, using the column-index.
  // That is: xcoord[nyNode][nxNode]
  // ... for efficiency (cpp row-major), loop over [nxNode] first.
  as3double (*xcoord)[nxNode] = (as3double (*)[nxNode]) GridCoordinate[0];
  as3double (*ycoord)[nxNode] = (as3double (*)[nxNode]) GridCoordinate[1];

  // Generate grid.
  as3double x = x0, y = y0;
  for(unsigned long j=0; j<nyNode; j++){
    // Reset the x-coordinate.
    x = x0;
    for(unsigned long i=0; i<nxNode; i++){

      // Book-keep the current grid coordinate.
      xcoord[j][i] = x;
      ycoord[j][i] = y;

      // Update the x-coordinate.
      x = x + dx;
    }
    // Update the y-coordinate.
    y = y + dy;
  }


  // // DEBUGGING
  // unsigned long idx = 0;
  // for(auto j=0; j<nyNode; j++){
  //   for(auto i=0; i<nxNode; i++){
  //
  //     unsigned long ij = j*nyNode+i;
  //
  //
  //     std::cout << xcoord[j][i] << ", "
  //               << ycoord[j][i] << ", "
  //               << GridCoordinate[0][idx] << ", "
  //               << GridCoordinate[1][idx] << ", "
  //               << GridCoordinate[0][ij] << ", "
  //               << GridCoordinate[1][ij] << std::endl;
  //     idx++;
  //   }
  // }

}



