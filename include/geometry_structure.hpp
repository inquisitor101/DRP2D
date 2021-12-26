#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"


class CGeometry {

  public:
    // Constructor.
    CGeometry(CConfig *config_container);

    // Destructor.
    ~CGeometry(void);

    // Getter: returns nPointCellP1.
		unsigned long GetnPointCellP1(void)                     const {return nPointCellP1;}
    // Getter: returns nxNode.
    unsigned long GetnxNode(void)                           const {return nxNode;}
    // Getter: returns nyNode.
    unsigned long GetnyNode(void)                           const {return nyNode;}
    // Getter: returns nNode.
    unsigned long GetnNode(void)                            const {return nNode;}
    // Getter: returns GridResolution[iDim].
    as3double GetGridResolution(unsigned short iDim)        const {return GridResolution[iDim];}
    // Getter: returns GridCoordinate[iDim].
    const as3double *GetGridCoordinate(unsigned short iDim) const {return GridCoordinate[iDim];}


	protected:

  private:
		// Total number of points in grid, for nPoly=1 cells.
		unsigned long nPointCellP1;

    // Number of nodes in x-direction.
    unsigned long nxNode;
    // Number of nodes in y-direction.
    unsigned long nyNode;
    // Total number of nodes.
    unsigned long nNode;

    // Grid domain size.
    // Convention: (WEST, EAST, SOUTH, NORTH).
    as3vector1d<as3double> GridSize;

    // Grid resolution.
    as3double GridResolution[nDim];

    // Grid coordinats at points in grid.
    // Dimension: [iDim][iNode].
    as3data1d<as3double> GridCoordinate;

    // Function that computes the total number of grid points based on the
    // assumption of nPoly=1 sub-elements. This is used to
    // compute nPointCellP1 for the VTK output.
    void ComputePointCellP1(void);

    // Function that generates a grid.
    void GenerateGrid(unsigned long nxNode,
                      unsigned long nyNode,
                      unsigned long nNode);
};


