#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"



class CInitial {

	public:
		// Constructor.
		CInitial(CConfig   *config_container,
						 CGeometry *geometry_container);

		// Destructor.
		virtual ~CInitial(void);

		// Pure virtual function that returns the initial condition in primitive form.
		// Must be overriden by one of the derived classes.
		virtual as3vector1d<as3double> SetInitialConditionPrimitive(as3double x,
                                                                as3double y,
                                                                as3double t) = 0;

    // Getter: returns betaPML.
    as3double GetBetaPML(void)                             const {return betaPML;}
    // Getter: returns kappax.
    virtual as3double GetKappax(void)                      const {return 0.0;}
    // Getter: returns kappat.
    virtual as3double GetKappat(void)                      const {return 0.0;}
    // Pure virtual getter: returns aInf. Must be implemented in a derived class.
    virtual as3double GetReferenceSpeedOfSound(void)       const = 0;
    // Pure virtual getter: returns V0. Must be implemented in a derived class.
    virtual as3double GetVelocityTransverse(void)          const = 0;
    // Pure virtual getter: returns U0. Must be implemented in a derived class.
    virtual as3double GetVelocityNormal(void)              const = 0;
    // Pure virtual getter: returns Tinf. Must be implemented in a derived class.
    virtual as3double GetTinf(void)                        const = 0;
    // Pure virtual getter: returns pInf. Must be implemented in a derived class.
    virtual as3double GetPinf(void)                        const = 0;
    // Pure virtual getter: returns A0. Must be implemented in a derived class.
    virtual as3double GetA0(void)                          const = 0;
    // Pure virtual getter: returns Mach. Must be implemented in a derived class.
    virtual as3double GetMach(void)                        const = 0;
    // Pure virtual getter: returns t0. Must be implemented in a derived class.
    virtual as3double Gett0(void)                          const = 0;
    // Pure virtual getter: returns true/false. Must be implemented in a derived class.
    virtual bool GetDimensionalProblem(void)               const = 0;
    // Pure virtual getter: returns Fxpseudo. Must be implemented in a derived class.
    virtual const as3vector1d<as3double> GetFxpseudo(void) const = 0;

	protected:
    // Dispersion-relaxation correction parameter. Used in PML formulation.
    as3double betaPML;

	private:

};


class CGaussianInitial : public CInitial {

	public:
		// Constructor.
		CGaussianInitial(CConfig   *config_container,
						 				 CGeometry *geometry_container);

		// Destructor.
		~CGaussianInitial(void) final;

		// Function that returns a Gaussian initial condition in primitive form.
    as3vector1d<as3double> SetInitialConditionPrimitive(as3double x,
                                                        as3double y,
                                                        as3double t) final;


    // Getter: returns aInf.
    as3double GetReferenceSpeedOfSound(void)       const final {return aInf;}
    // Getter: returns V0.
    as3double GetVelocityTransverse(void)          const final {return vInf;}
    // Getter: returns U0.
    as3double GetVelocityNormal(void)              const final {return uInf;}
    // Getter: returns Tinf.
    as3double GetTinf(void)                        const final {return Tinf;}
    // Getter: returns pInf.
    as3double GetPinf(void)                        const final {return pInf;}
    // Getter: returns A0.
    as3double GetA0(void)                          const final {return A0;}
    // Getter: returns Mach.
    as3double GetMach(void)                        const final {return Mach;}
    // Getter: returns t0. Not needed in this IC.
    as3double Gett0(void)                          const final {return 0.0;}
    // Getter: returns true/false.
    bool GetDimensionalProblem(void)               const final {return true;}
    // Getter: returns Fxpseudo.
    const as3vector1d<as3double> GetFxpseudo(void) const final {return Fxpseudo;}

	protected:

	private:
		// Pulse center.
		as3double x0;
		as3double y0;
		// Pulse strength.
		as3double A0;
		// Pulse width.
		as3double b;
		// Pulse attenuation.
		as3double kappa;
    // Background Mach number.
    as3double Mach;
    // Flow direction [degrees].
    as3double theta;
    // Velocity magnitude free-stream.
    as3double VelInf;
		// Free-stream values.
		as3double uInf;
		as3double vInf;
    as3double pInf;
    as3double Tinf;
		as3double rhoInf;
    as3double aInf;
    // Flux in x-dimension based on pseudo-mean flow.
    as3vector1d<as3double> Fxpseudo;
};


class CIsentropicVortexInitial : public CInitial {

	public:
		// Constructor.
		CIsentropicVortexInitial(CConfig   *config_container,
        						 				 CGeometry *geometry_container);

		// Destructor.
		~CIsentropicVortexInitial(void) final;

		// Function that returns an isentropic vortex initial condition in primitive form.
    as3vector1d<as3double> SetInitialConditionPrimitive(as3double x,
                                                        as3double y,
                                                        as3double t) final;

    // Getter: returns aInf.
    as3double GetReferenceSpeedOfSound(void)       const final {return aInf;}
    // Getter: returns V0.
    as3double GetVelocityTransverse(void)          const final {return vInf;}
    // Getter: returns U0.
    as3double GetVelocityNormal(void)              const final {return uInf;}
    // Getter: returns Tinf.
    as3double GetTinf(void)                        const final {return Tinf;}
    // Getter: returns pInf.
    as3double GetPinf(void)                        const final {return pInf;}
    // Getter: returns A0.
    as3double GetA0(void)                          const final {return A0;}
    // Getter: returns Mach.
    as3double GetMach(void)                        const final {return Mach;}
    // Getter: returns t0. Not needed in this IC.
    as3double Gett0(void)                          const final {return 0.0;}
    // Getter: returns DimensionalProblem.
    bool GetDimensionalProblem(void)               const final {return true;}
    // Getter: returns Fxpseudo.
    const as3vector1d<as3double> GetFxpseudo(void) const final {return Fxpseudo;}

	protected:

	private:
		// Vortex center.
		as3double x0;
		as3double y0;
		// Vortex strength.
		as3double A0;
		// Vortex radius.
		as3double Rv;
    // Background Mach number.
    as3double Mach;
    // Flow direction [degrees].
    as3double theta;
    // Velocity magnitude free-stream.
    as3double VelInf;
    // Abbreviations.
    as3double ovRv2;
    as3double alpha;
    as3double omega;
		// Free-stream values.
		as3double uInf;
		as3double vInf;
    as3double pInf;
    as3double Tinf;
		as3double rhoInf;
    as3double aInf;
    // Flux in x-dimension based on pseudo-mean flow.
    as3vector1d<as3double> Fxpseudo;
};
