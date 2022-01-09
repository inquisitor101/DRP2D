#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "initial_structure.hpp"
#include "stencil_structure.hpp"



class CSpatial {

	public:
		// Constructor.
		CSpatial(CConfig   	   *config_container,
						 CGeometry 	   *geometry_container,
             CInitial      *initial_container,
             CStencil     **stencil_container,
             unsigned short iZone);

		// Destructor.
		virtual ~CSpatial(void);

    // Pure virtual function that sweeps through the grid and updates the
    // residual accordingly. Must be overriden by a derived class.
    virtual void ComputeResidual(CConfig                *config_container,
                                 CGeometry              *geometry_container,
                                 CSpatial               *spatial_container,
                                 CInitial               *initial_container,
                                 as3data1d<as3double>   &work_array,
                                 as3data1d<as3double>   &residual,
                                 as3double               localTime,
                                 as3vector1d<as3double> &MonitoringData) = 0;

  protected:
    // Current zone ID.
    unsigned short zoneID;

    // Number of nodes in x-direction.
    unsigned long nxNode;
    // Number of nodes in y-direction.
    unsigned long nyNode;
		// Total number of nodes.
		unsigned long nNode;

    // Total nodes in x-direction, including ghost nodes.
    unsigned long nxTotal;
    // Total nodes in y-direction, including ghost nodes.
    unsigned long nyTotal;

    // Stencil size in west(left) x-direction.
		unsigned short NxStencil;
		// Stencil size in east(right) x-direction.
		unsigned short MxStencil;

		// Stencil size in south(bottom) y-direction.
		unsigned short NyStencil;
		// Stencil size in north(top) y-direction.
		unsigned short MyStencil;

    // Non-zero stencil indices in both dimensions.
    as3vector2d<long> ijStencil;
    // Non-zero stencil coefficients in both dimensions.
    as3vector2d<as3double> cStencil;

    // Indices of all nodes needed, included ghost ones. This excludes the
    // corner ghost nodes that are unecessary annd might result in floating
    // point issues, such as inf and nan.
    as3vector1d<unsigned long> RequiredIndices;

  private:
    // Function that initializes the finite-difference stencil.
    void InitializeStencilStrategy(CConfig   *config_container,
                                   CGeometry *geometry_container);

    // Function that identifies the required nodal indices.
    void IdentifyRequiredNodalIndices(CConfig   *config_container,
                                      CGeometry *geometry_container);

};


class CEESpatial : public CSpatial {

	public:
		// Constructor.
		CEESpatial(CConfig   	   *config_container,
  						 CGeometry 	   *geometry_container,
               CInitial      *initial_container,
               CStencil     **stencil_container,
               unsigned short iZone);

		// Destructor.
		~CEESpatial(void);

  protected:
    // Function that sweeps through the grid and updates the residual accordingly.
    void ComputeResidual(CConfig                *config_container,
                         CGeometry              *geometry_container,
                         CSpatial               *spatial_container,
                         CInitial               *initial_container,
                         as3data1d<as3double>   &work_array,
                         as3data1d<as3double>   &residual,
                         as3double               localTime,
                         as3vector1d<as3double> &MonitoringData);

  private:
};



