#include "geometry_structure.hpp"



CGeometry::CGeometry
(
 CConfig *config_container
)
 /*
  * Constructor, Reads input grid.
  */
{
  // Set number of zones.
  nZone = config_container->GetnZone();

	// Reserve memory for geometry zone container.
	geometry_zone.resize(nZone, nullptr);

  // Create each zone.
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    geometry_zone[iZone] = new CGeometryZone(config_container, iZone);

  // Compute the total number of nPointCellP1.
  ComputePointCellP1();
}


CGeometry::~CGeometry
(
 void
)
 /*
  * Destructor for CGeometry class, frees allocated memory.
  */
{
  for(unsigned short i=0; i<geometry_zone.size(); i++)
    if( geometry_zone[i] ) delete geometry_zone[i];
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
  nPointCellP1 = 0;
  for(unsigned short iZone=0; iZone<nZone; iZone++)
    nPointCellP1 += geometry_zone[iZone]->GetnCell();

  // Include the quadrilateral vertices.
  nPointCellP1 *= N_POINTS_QUADRILATERAL;
}


CGeometryZone::CGeometryZone
(
 CConfig       *config_container,
 unsigned short iZone
)
 /*
  * Constructor, Reads input grid.
  */
{
  // Set current zone ID.
  zoneID = iZone;

  // Extract number of nodes in x-direction.
  nxNode = config_container->GetnxNode(iZone);
  // Extract number of nodes in y-direction.
  nyNode = config_container->GetnyNode(iZone);
	// Deduce total number of nodes.
	nNode  = nxNode*nyNode;

  // Number of cells in x-direction.
  nxCell = nxNode-1;
  // Number of cells in y-direction.
  nyCell = nyNode-1;
  // Deduce total number of cells.
  nCell  = nxCell*nyCell;

  // Generate grid zone.
  GenerateGrid(config_container);
}


CGeometryZone::~CGeometryZone
(
 void
)
 /*
  * Destructor for CGeometryZone class, frees allocated memory.
  */
{
  for(unsigned short i=0; i<GridCoordinate.size(); i++)
    if( GridCoordinate[i] ) delete [] GridCoordinate[i];
}



void CGeometryZone::GenerateGrid
(
 CConfig *config_container
)
 /*
  * Function that generates a grid in current zone.
  */
{
  // Main/physical domain bounding box.
  auto MainBox = config_container->GetDomainBound();

  // Partition the main domain explicitly for readability.
  const as3double xmin0 = MainBox[0];
  const as3double xmax0 = MainBox[1];
  const as3double ymin0 = MainBox[2];
  const as3double ymax0 = MainBox[3];

  // Number of nodes in x-direction in main zone.
  const unsigned long nx0 = config_container->GetnxNode(ZONE_MAIN);
  // Number of nodes in y-direction in main zone.
  const unsigned long ny0 = config_container->GetnyNode(ZONE_MAIN);

  // Domain size in x-direction of main zone.
  const as3double lx0 = xmax0 - xmin0;
  // Domain size in y-direction of main zone.
  const as3double ly0 = ymax0 - ymin0;

  // Element size in x-direction in the main zone.
  const as3double dx0 = lx0 / (as3double) nx0;
  // Element size in y-direction in the main zone.
  const as3double dy0 = ly0 / (as3double) ny0;

  // Number of cells in x-direction in main zone.
  const unsigned long nxCell0 = nx0-1;
  // Number of cells in y-direction in main zone.
  const unsigned long nyCell0 = ny0-1;

  // Temporary storage for the main grid resolution.
  as3vector1d<as3double> hx0(nxCell0);
  as3vector1d<as3double> hy0(nyCell0);

  // Check whether this is a uniform grid or not.
  if( config_container->GetUniformGridResolution() ){

    // Compute step-size in x-direction.
    for(unsigned long i=0; i<nxCell0; i++) hx0[i] = dx0;
    // Compute step-size in y-direction.
    for(unsigned long j=0; j<nyCell0; j++) hy0[j] = dy0;
  }
  else {
    // This is not a uniform grid resolution, not support yet.
    Terminate("CGeometryZone::GenerateGrid", __FILE__, __LINE__,
              "Only uniform grids are supported for now.");
  }

  // Allocate temporary storage for generic grid resolution.
  as3vector1d<as3double> hx(nxCell);
  as3vector1d<as3double> hy(nyCell);

  // Determine what zone we are dealing with.
  switch( config_container->GetTypeZone(zoneID) ){

    // Main zone.
    case(ZONE_MAIN):
    {
      // Copy the pre-computed main zone resolution.
      hx = hx0; hy = hy0;

      // Determine domain size.
      ZoneSize = MainBox;

      break;
    }

    // West zone.
    case(ZONE_WEST):
    {
      // Main zone cell indices.
      unsigned long II = 0;

      // Expansion ratio in x-dimension in west zone.
      const as3double rx = config_container->GetdsNodalRatioZone()[0];

      // Compute step-size in x-direction.
      for(unsigned long i=0; i<nxCell; i++) hx[i] = hx0[II]*pow(rx, nxNode-i);
      // Copy the pre-compute step-size in y-direction in main zone.
      hy = hy0;

      // Offset between main zone and west zone in x-direction.
      const as3double dhx = hx0[II]*rx;
      // Total zone width in x-dimension.
      as3double ww = dhx; for(auto ds : hx) ww += ds;

      // Determine domain size.
      ZoneSize[0] = xmin0 - ww;  // xmin.
      ZoneSize[1] = xmin0 - dhx; // xmax.
      ZoneSize[2] = ymin0;       // ymin.
      ZoneSize[3] = ymax0;       // ymax.

      break;
    }

    // East zone.
    case(ZONE_EAST):
    {
      // Main zone cell indices.
      unsigned long II = nyCell0-1;

      // Expansion ratio in x-dimension in east zone.
      const as3double rx = config_container->GetdsNodalRatioZone()[1];

      // Compute step-size in x-direction.
      for(unsigned long i=0; i<nxCell; i++) hx[i] = hx0[II]*pow(rx, i+2);
      // Copy the pre-compute step-size in y-direction in main zone.
      hy = hy0;

      // Offset between main zone and east zone in x-direction.
      const as3double dhx = hx0[II]*rx;
      // Total zone width in x-dimension.
      as3double ww = dhx; for(auto ds : hx) ww += ds;

      // Determine domain size.
      ZoneSize[0] = xmax0 + dhx; // xmin.
      ZoneSize[1] = xmax0 + ww;  // xmax.
      ZoneSize[2] = ymin0;       // ymin.
      ZoneSize[3] = ymax0;       // ymax.

      break;
    }

    // South zone.
    case(ZONE_SOUTH):
    {
      // Main zone cell indices.
      unsigned long JJ = 0;

      // Expansion ratio in y-dimension in south zone.
      const as3double ry = config_container->GetdsNodalRatioZone()[2];

      // Copy the pre-compute step-size in x-direction in main zone.
      hx = hx0;
      // Compute step-size in y-direction.
      for(unsigned long j=0; j<nyCell; j++) hy[j] = hy0[JJ]*pow(ry, nyNode-j);

      // Offset between main zone and south zone in y-direction.
      const as3double dhy = hy0[JJ]*ry;
      // Total zone width in y-dimension.
      as3double ww = dhy; for(auto ds : hy) ww += ds;

      // Determine domain size.
      ZoneSize[0] = xmin0;       // xmin.
      ZoneSize[1] = xmax0;       // xmax.
      ZoneSize[2] = ymin0 - ww;  // ymin.
      ZoneSize[3] = ymin0 - dhy; // ymax.

      break;
    }

    // North zone.
    case(ZONE_NORTH):
    {
      // Main zone cell indices.
      unsigned long JJ = nyCell0-1;

      // Expansion ratio in y-dimension in north zone.
      const as3double ry = config_container->GetdsNodalRatioZone()[3];

      // Copy the pre-compute step-size in x-direction in main zone.
      hx = hx0;
      // Compute step-size in y-direction.
      for(unsigned long j=0; j<nyCell; j++) hy[j] = hy0[JJ]*pow(ry, j+2);

      // Offset between main zone and north zone in x-direction.
      const as3double dhy = hy0[JJ]*ry;
      // Total zone width in y-dimension.
      as3double ww = dhy; for(auto ds : hy) ww += ds;

      // Determine domain size.
      ZoneSize[0] = xmin0;       // xmin.
      ZoneSize[1] = xmax0;       // xmax.
      ZoneSize[2] = ymax0 + dhy; // ymin.
      ZoneSize[3] = ymax0 + ww;  // ymax.

      break;
    }

    // Corner0 zone.
    case(ZONE_CORNER_0):
    {
      // Main zone cell indices.
      unsigned long II = 0, JJ = 0;

      // Expansion ratio in x-dimension in west zone.
      const as3double rx = config_container->GetdsNodalRatioZone()[0];
      // Expansion ratio in y-dimension in south zone.
      const as3double ry = config_container->GetdsNodalRatioZone()[2];

      // Compute step-size in x-direction.
      for(unsigned long i=0; i<nxCell; i++) hx[i] = hx0[II]*pow(rx, nxNode-i);
      // Compute step-size in y-direction.
      for(unsigned long j=0; j<nyCell; j++) hy[j] = hy0[JJ]*pow(ry, nyNode-j);

      // Offset between main zone and west zone in x-direction.
      const as3double dhx = hx0[II]*rx;
      // Offset between main zone and south zone in y-direction.
      const as3double dhy = hy0[JJ]*ry;

      // Total zone width in x-dimension.
      as3double wwx = dhx; for(auto ds : hx) wwx += ds;
      // Total zone width in y-dimension.
      as3double wwy = dhy; for(auto ds : hy) wwy += ds;

      // Determine domain size.
      ZoneSize[0] = xmin0 - wwx; // xmin.
      ZoneSize[1] = xmin0 - dhx; // xmax.
      ZoneSize[2] = ymin0 - wwy; // ymin.
      ZoneSize[3] = ymin0 - dhy; // ymax.

      break;
    }

    // Corner1 zone.
    case(ZONE_CORNER_1):
    {
      // Main zone cell indices.
      unsigned long II = nxCell0-1, JJ = 0;

      // Expansion ratio in x-dimension in east zone.
      const as3double rx = config_container->GetdsNodalRatioZone()[1];
      // Expansion ratio in y-dimension in south zone.
      const as3double ry = config_container->GetdsNodalRatioZone()[2];

      // Compute step-size in x-direction.
      for(unsigned long i=0; i<nxCell; i++) hx[i] = hx0[II]*pow(rx, i+2);
      // Compute step-size in y-direction.
      for(unsigned long j=0; j<nyCell; j++) hy[j] = hy0[JJ]*pow(ry, nyNode-j);

      // Offset between main zone and east zone in x-direction.
      const as3double dhx = hx0[II]*rx;
      // Offset between main zone and south zone in y-direction.
      const as3double dhy = hy0[JJ]*ry;

      // Total zone width in x-dimension.
      as3double wwx = dhx; for(auto ds : hx) wwx += ds;
      // Total zone width in y-dimension.
      as3double wwy = dhy; for(auto ds : hy) wwy += ds;

      // Determine domain size.
      ZoneSize[0] = xmax0 + dhx; // xmin.
      ZoneSize[1] = xmax0 + wwx; // xmax.
      ZoneSize[2] = ymin0 - wwy; // ymin.
      ZoneSize[3] = ymin0 - dhy; // ymax.

      break;
    }

    // Corner2 zone.
    case(ZONE_CORNER_2):
    {
      // Main zone cell indices.
      unsigned long II = 0, JJ = nyCell0-1;

      // Expansion ratio in x-dimension in west zone.
      const as3double rx = config_container->GetdsNodalRatioZone()[0];
      // Expansion ratio in y-dimension in north zone.
      const as3double ry = config_container->GetdsNodalRatioZone()[3];

      // Compute step-size in x-direction.
      for(unsigned long i=0; i<nxCell; i++) hx[i] = hx0[II]*pow(rx, nxNode-i);
      // Compute step-size in y-direction.
      for(unsigned long j=0; j<nyCell; j++) hy[j] = hy0[JJ]*pow(ry, j+2);

      // Offset between main zone and west zone in x-direction.
      const as3double dhx = hx0[II]*rx;
      // Offset between main zone and north zone in y-direction.
      const as3double dhy = hy0[JJ]*ry;

      // Total zone width in x-dimension.
      as3double wwx = dhx; for(auto ds : hx) wwx += ds;
      // Total zone width in y-dimension.
      as3double wwy = dhy; for(auto ds : hy) wwy += ds;

      // Determine domain size.
      ZoneSize[0] = xmin0 - wwx; // xmin.
      ZoneSize[1] = xmin0 - dhx; // xmax.
      ZoneSize[2] = ymax0 + dhy; // ymin.
      ZoneSize[3] = ymax0 + wwy; // ymax.

      break;
    }

    // Corner3 zone.
    case(ZONE_CORNER_3):
    {
      // Main zone cell indices.
      unsigned long II = nxCell0-1, JJ = nyCell0-1;

      // Expansion ratio in x-dimension in east zone.
      const as3double rx = config_container->GetdsNodalRatioZone()[1];
      // Expansion ratio in y-dimension in north zone.
      const as3double ry = config_container->GetdsNodalRatioZone()[3];

      // Compute step-size in x-direction.
      for(unsigned long i=0; i<nxCell; i++) hx[i] = hx0[II]*pow(rx, i+2);
      // Compute step-size in y-direction.
      for(unsigned long j=0; j<nyCell; j++) hy[j] = hy0[JJ]*pow(ry, j+2);

      // Offset between main zone and east zone in x-direction.
      const as3double dhx = hx0[II]*rx;
      // Offset between main zone and north zone in y-direction.
      const as3double dhy = hy0[JJ]*ry;

      // Total zone width in x-dimension.
      as3double wwx = dhx; for(auto ds : hx) wwx += ds;
      // Total zone width in y-dimension.
      as3double wwy = dhy; for(auto ds : hy) wwy += ds;

      // Determine domain size.
      ZoneSize[0] = xmax0 + dhx; // xmin.
      ZoneSize[1] = xmax0 + wwx; // xmax.
      ZoneSize[2] = ymax0 + dhy; // ymin.
      ZoneSize[3] = ymax0 + wwy; // ymax.

      break;
    }

    default:
      Terminate("CGeometryZone::GenerateGridZone", __FILE__, __LINE__,
                "Unknown/wrong zone type detected");
  }


  // Cast 1D array into 2D sub-arrays for readability, using the column-index.
  // That is: xcoord[nyNode][nxNode]
  // ... for efficiency (cpp row-major), loop over [nxNode] first.
  as3double (*xcoord)[nxNode] = (as3double (*)[nxNode]) GridCoordinate[0];
  as3double (*ycoord)[nxNode] = (as3double (*)[nxNode]) GridCoordinate[1];

  // Initialize ymin grid nodes.
  for(unsigned long i=0; i<nxNode; i++) ycoord[i][0] = ZoneSize[2]; // y = ymin
  // Initialize xmin grid nodes.
  for(unsigned long j=0; j<nyNode; j++) xcoord[0][j] = ZoneSize[0]; // x = xmin

  // Generate grid.
  for(unsigned long j=0; j<nyCell; j++){
    for(unsigned long i=0; i<nxCell; i++){

      // Book-keep the current grid coordinate.
      xcoord[j][i+1] = xcoord[j][i] + hx[i];
      ycoord[j+1][i] = ycoord[j][i] + hy[j];
    }
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
