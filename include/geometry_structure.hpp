#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometryZone;


class CGeometry {

  public:
    // Constructor.
    CGeometry(CConfig *config_container);

    // Destructor.
    ~CGeometry(void);

    // Getter: returns nPointCellP1.
		unsigned long GetnPointCellP1(void) const {return nPointCellP1;}

	protected:

  private:
    // Number of zones.
    unsigned short nZone;
		// Total number of points in grid, for nPoly=1 cells.
		unsigned long nPointCellP1;

    // Grid geometry in each zone.
		as3data1d<CGeometryZone> geometry_zone;

    // Function that computes the total number of grid points based on the
    // assumption of nPoly=1 sub-elements. This is used to
    // compute nPointCellP1 for the VTK output.
    void ComputePointCellP1(void);
};



class CGeometryZone {

  public:
    // Constructor.
    CGeometryZone(CConfig       *config_container,
                  unsigned short iZone);

    // Destructor.
    ~CGeometryZone(void);

    // Getter: returns nxNode.
    unsigned long GetnxNode(void)                           const {return nxNode;}
    // Getter: returns nyNode.
    unsigned long GetnyNode(void)                           const {return nyNode;}
    // Getter: returns nNode.
    unsigned long GetnNode(void)                            const {return nNode;}
    // Getter: returns nxCell.
    unsigned long GetnxCell(void)                           const {return nxCell;}
    // Getter: returns nyCell.
    unsigned long GetnyCell(void)                           const {return nyCell;}
    // Getter: returns nCell.
    unsigned long GetnCell(void)                            const {return nCell;}
    // Getter: returns GridResolution[iDim].
    as3double GetGridResolution(unsigned short iDim)        const {return GridResolution[iDim];}
    // Getter: returns GridCoordinate[iDim].
    const as3double *GetGridCoordinate(unsigned short iDim) const {return GridCoordinate[iDim];}

  protected:

  private:
    // Current zone ID.
    unsigned short zoneID;
    // Number of nodes in x-direction.
    unsigned long nxNode;
    // Number of nodes in y-direction.
    unsigned long nyNode;
    // Total number of nodes.
    unsigned long nNode;

    // Number of cells in x-direction.
    unsigned long nxCell;
    // Number of cells in y-direction.
    unsigned long nyCell;
    // Total number of cells.
    unsigned long nCell;

    // Grid domain size in this zone.
    // Convention: (WEST, EAST, SOUTH, NORTH).
    as3vector1d<as3double> ZoneSize;

    // Grid resolution.
    as3double GridResolution[nDim];

    // Grid coordinats at points in grid.
    // Dimension: [iDim][iNode].
    as3data1d<as3double> GridCoordinate;

    // Function that generates a grid.
    void GenerateGrid(CConfig *config_container);
};
