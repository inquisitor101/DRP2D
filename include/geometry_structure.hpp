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
		unsigned long GetnPointCellP1(void)                        const {return nPointCellP1;}
    // Getter: returns nxNode.
    unsigned long GetnxNode(void)                              const {return nxNode;}
    // Getter: returns nyNode.
    unsigned long GetnyNode(void)                              const {return nyNode;}
    // Getter: returns nNode.
    unsigned long GetnNode(void)                               const {return nNode;}
    // Getter: returns nNodeBuffer.
    unsigned long GetnNodeBuffer(void)                         const {return nNodeBuffer;}
    // Getter: returns GridResolution[iDim].
    as3double GetGridResolution(unsigned short iDim)           const {return GridResolution[iDim];}
    // Getter: returns GridCoordinate[iDim].
    const as3double *GetGridCoordinate(unsigned short iDim)    const {return GridCoordinate[iDim];}
    // Getter: returns IBufferLayer.
    const as3vector1d<unsigned long> &GetIBufferLayer(void)    const {return IBufferLayer;}
    // Getter: returns DampingFunctionXDir.
    const as3vector1d<as3double> &GetDampingFunctionXDir(void) const {return DampingFunctionXDir;}

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

    // Number of nodes in buffer layer in x-direction.
    unsigned long nxNodeBuffer;
    // Number of nodes in buffer layer in y-direction.
    unsigned long nyNodeBuffer;
    // Total number of nodes in buffer layer.
    unsigned long nNodeBuffer;

    // Buffer layer i-indices.
    as3vector1d<unsigned long> IBufferLayer;
    // Buffer layer damping in x-direction.
    as3vector1d<as3double> DampingFunctionXDir;

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

    // Function that generates the buffer layer.
    void GenerateBufferLayer(CConfig *config_container);
};


