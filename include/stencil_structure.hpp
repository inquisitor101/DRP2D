#pragma once

#include "option_structure.hpp"
#include "geometry_structure.hpp"
#include "config_structure.hpp"



class CStencil {

	public:
		// Constructor.
		CStencil(CConfig   *config_container,
					   CGeometry *geometry_container);

		// Destructor.
		virtual ~CStencil(void);

    // Pure virtual getter: returns TypeStencil. Must be implemented in a derived class.
    virtual unsigned short GetTypeStencil(void)               const = 0;
    // Pure virtual getter: returns cStencil.size(). Must be implemented in a derived class.
    virtual unsigned short GetnStencil(void)                  const = 0;
    // Pure virtual getter: returns MMStencil. Must be implemented in a derived class.
    virtual unsigned short GetMMStencil(void)                 const = 0;
    // Pure virtual getter: returns NNStencil. Must be implemented in a derived class.
    virtual unsigned short GetNNStencil(void)                 const = 0;
    // Pure virtual getter: returns kStencil. Must be implemented in a derived class.
    virtual const as3vector1d<short>     GetkStencil(void)    const = 0;
    // Pure virtual getter: returns cStencil. Must be implemented in a derived class.
    virtual const as3vector1d<as3double> GetcStencil(void)    const = 0;

	protected:

	private:

};


class CDRPM3N3Stencil : public CStencil {

	public:
		// Constructor.
		CDRPM3N3Stencil(CConfig   *config_container,
      					    CGeometry *geometry_container);

		// Destructor.
		~CDRPM3N3Stencil(void) final;

    // Getter: returns TypeStencil.
    unsigned short GetTypeStencil(void)               const {return STENCIL_DRP_M3N3;}
    // Getter: returns cStencil.size().
    unsigned short GetnStencil(void)                  const {return cStencil.size();}
    // Getter: returns MMStencil.
    unsigned short GetMMStencil(void)                 const {return MMStencil;}
    // Getter: returns NNStencil.
    unsigned short GetNNStencil(void)                 const {return NNStencil;}
    // Getter: returns kStencil.
    const as3vector1d<short>     GetkStencil(void)    const {return kStencil;}
    // Getter: returns cStencil.
    const as3vector1d<as3double> GetcStencil(void)    const {return cStencil;}

	protected:

	private:
    // Stencil size in backward direction.
    unsigned short NNStencil;
    // Stencil size in forward direction.
    unsigned short MMStencil;

    // Stencil indices offset.
    as3vector1d<short> kStencil;
    // Stencil coefficients.
    as3vector1d<as3double> cStencil;
};





