#pragma once

#include "option_structure.hpp"
#include <map>



class CConfig {

	public:
		// Constructor.
		CConfig(const char *configFile);

		// Destructor.
		~CConfig(void);

    // Getter: returns nxNode.
    unsigned long GetnxNode(void)                              const {return nxNode;}
    // Getter: returns nyNode.
    unsigned long GetnyNode(void)                              const {return nyNode;}
    // Getter: returns NStencil[iDim].
    unsigned short GetnnStencil(unsigned short iDim)           const {return NStencil[iDim];}
    // Getter: returns MStencil[iDim].
    unsigned short GetmmStencil(unsigned short iDim)           const {return MStencil[iDim];}
    // Getter: returns NxStencil.
    unsigned short GetnxStencil(void)                          const {return NStencil[0];}
    // Getter: returns MxStencil.
    unsigned short GetmxStencil(void)                          const {return MStencil[0];}
    // Getter: returns NyStencil.
    unsigned short GetnyStencil(void)                          const {return NStencil[1];}
    // Getter: returns MyStencil.
    unsigned short GetmyStencil(void)                          const {return MStencil[1];}
    // Getter: returns SimulationTime[0].
		as3double GetSimulationStartTime(void)                     const {return SimulationTime[0];}
		// Getter: returns SimulationTime[1].
		as3double GetSimulationFinalTime(void)                     const {return SimulationTime[1];}
    // Getter: returns MaxIter.
    unsigned long GetMaxIter(void)                             const {return MaxIter;}
    // Getter: returns OutputVTKFilename.
    const char *GetOutputVTKFilename(void)                     const {return OutputVTKFilename.c_str();}
    // Getter: returns WriteFreq.
    unsigned long GetWriteFreq(void)                           const {return WriteFreq;}
    // Getter: returns OutputFreq.
    unsigned long GetOutputFreq(void)                          const {return OutputFreq;}
    // Getter: returns DomainBound.
    const as3vector1d<as3double> GetDomainBound(void)          const {return DomainBound;}
    // Getter: returns MachInf.
    as3double GetMachInf(void)                                 const {return MachInf;}
    // Getter: returns FlowAngle.
    as3double GetFlowAngle(void)                               const {return FlowAngle;}
    // Getter: returns CrossFlow.
    bool GetCrossFlow(void)                                    const {return CrossFlow;}
    // Getter: returns CenterX0.
    const as3vector1d<as3double> &GetCenterX0(void)            const {return CenterX0;}
    // Getter: returns DisturbanceRatio.
    as3double GetDisturbanceRatio(void)                        const {return DisturbanceRatio;}
    // Getter: returns DisturbanceWidth.
    as3double GetDisturbanceWidth(void)                        const {return DisturbanceWidth;}
    // Getter: returns TypeSolver.
    unsigned short GetTypeSolver(void)                         const {return TypeSolver;}
    // Getter: returns TypeIC.
    unsigned short GetTypeIC(void)                             const {return TypeIC;}
    // Getter: returns TypeTemporalScheme.
    unsigned short GetTypeTemporalScheme(void)                 const {return TypeTemporalScheme;}
    // Getter: returns TypeExternalBC[iBoundary].
    unsigned short GetTypeExternalBC(unsigned short iBoundary) const {return TypeExternalBC[iBoundary];}
    // Getter: returns CFL.
    as3double GetCFL(void)                                     const {return CFL;}
    // Getter: returns TimeStep.
    as3double GetTimeStep(void)                                const {return TimeStep;}
    // Getter: returns AdaptTime.
    bool GetAdaptTime(void)                                    const {return AdaptTime;}
    // Getter: returns TypeStencil[iDim].
    unsigned short GetTypeStencil(unsigned short iDim)         const {return TypeStencil[iDim];}
    // Getter: returns NodeBufferLayer.
    const as3vector1d<unsigned long> &GetNodeBufferLayer(void) const {return NodeBufferLayer;}

	private:
    // Default data parameters.
		struct {
			// Restart solution or not.
			std::string NameRestartSolution = "false";
			// Time step.
			as3double TimeStep = -1.0;
			// Maximum iterations in time.
			unsigned long MaxIter = 100000;
      // Output writing frequency.
      unsigned long WriteFreq  = 10;
      // Output monitoring frequency.
      unsigned long OutputFreq = 100;
      // CFL number.
      as3double CFL = 1.0;
      // Adaptive time step.
      std::string NameAdaptTime = "false";
      // Type of buffer layer.
      as3vector1d<std::string> NameTypeBufferLayer = { "NONE" };
      // Sponge-layer damping exponential coefficient.
      as3vector1d<as3double> DampingExponent = { 2.0 };
      // Sponge-layer damping constant.
      as3vector1d<as3double> DampingConstant = { 2.0 };
      // Grid-stretching constant.
      as3vector1d<as3double> GridStretchingConstant = { 1.0 };
      // Grid-stretching exponent.
      as3vector1d<as3double> GridStretchingExponent = { 1.0 };
      // Number of nodes per buffer layer.
      as3vector1d<unsigned long> NodeBufferLayer = { 10 };
      // Free-stream Mach number.
      as3double MachInf = 0.5;
      // Flow angle in degrees (w.r.t. x-direction).
      as3double FlowAngle = 0.0;
      // Disturbance center.
      as3vector1d<as3double> CenterX0 = { 0.0, 0.0 };
      // Ratio of disturbance w.r.t. background flow.
      as3double DisturbanceRatio = 0.2;
      // Disturbance width.
      as3double DisturbanceWidth = 0.25;
      // Grid-stretching used.
      bool GridStretching = false;

		} DefaultParam;


		// Number of physical nodes in x-direction.
		unsigned long nxNode;
		// Number of physical nodes in y-direction.
		unsigned long nyNode;

    // Stencil size in backward-direction, w.r.t. grid orientation.
    unsigned short NStencil[nDim];
    // Stencil size in forward-direction, w.r.t. grid orientation.
    unsigned short MStencil[nDim];

    // Type of stencils used in eaech direction, string name.
    as3vector1d<std::string> NameTypeStencil;
    // Type of stencil used in x-direction.
    as3vector1d<unsigned short> TypeStencil;

    // Domain bounds (physical domain).
    // Convention: (WEST, EAST, SOUTH, NORTH).
    as3vector1d<as3double> DomainBound;
    // Type of solver per zone.
    unsigned short TypeSolver;

    // Type of buffer layer, string name.
    as3vector1d<std::string> NameTypeBufferLayer;
    // Type of buffer layer, enum values.
    as3vector1d<unsigned short> TypeBufferLayer;

    // Decide whether or not to restart a simulation.
    bool RestartSolution;
    // Output writing frequency.
    unsigned long WriteFreq;
    // Output monitoring frequency.
    unsigned long OutputFreq;
    // Output solution filename.
    std::string OutputSolFilename;
    // Output visualization filename.
    std::string OutputVTKFilename;
    // Restart solution filename.
    std::string RestartFilename;

    // Type of temporal discretization, string name.
    std::string NameTemporalScheme;
    // Type of temporal discretization.
    unsigned short TypeTemporalScheme;
    // Simulation start and end times.
    // index: [0]: start, [1]: end.
    as3double SimulationTime[2];
    // Time step input.
    as3double TimeStep;
    // Maximum temporal iterations.
    unsigned long MaxIter;
    // CFL number.
    as3double CFL;
    // Adaptive time step.
    bool AdaptTime;

    // Free-stream Mach number.
    as3double MachInf;
    // Flow angle in degrees (w.r.t. x-direction).
    as3double FlowAngle;
    // Presence of a cross-flow.
    bool CrossFlow;
    // Disturbance center.
    as3vector1d<as3double> CenterX0;
    // Ratio of disturbance w.r.t. background flow.
    as3double DisturbanceRatio;
    // Disturbance width.
    as3double DisturbanceWidth;

    // Type of initial conditions, string name.
    std::string NameInitialCondition;
    // Type of initial conditions, enum values.
    unsigned short TypeIC;

    // Sponge-layer damping exponential coefficient.
    as3vector1d<as3double> DampingExponent;
    // Sponge-layer damping constant.
    as3vector1d<as3double> DampingConstant;
    // Number of nodes per buffer layer.
    as3vector1d<unsigned long> NodeBufferLayer;

    // Whether or not to use grid-stretching.
    bool GridStretching;
    // Grid-stretching constant.
    as3vector1d<as3double> GridStretchingConstant;
    // Grid-stretching exponent.
    as3vector1d<as3double> GridStretchingExponent;

    // Type of boundary conditions, string name:
    // ... indices: (SOUTH, NORTH, WEST, EAST).
    as3vector1d<std::string> NameBoundaryCondition;
    // Type of external boundary conditions over all zones:
    // ... indices: (SOUTH, NORTH, WEST, EAST).
    as3vector1d<unsigned short> TypeExternalBC;


    // Function that maps TemporalScheme from string to enum.
    void MapTemporalScheme(void);
    // Function that maps TypeIC from string to enum.
    void MapTypeIC(void);
    // Function that maps TypeExternalBC from string to enum.
    void MapTypeExternalBC(void);
    // Function that maps TypeStencil from string to enum.
    void MapTypeStencil(void);
    // Function that maps TypeBufferLayer from string to enum.
    void MapTypeBufferLayer(void);


    // Function that reads grid options.
    bool ReadGridOptions(const char *configFile);
    // Function that reads input/output options.
  	bool ReadIOOptions(const char *configFile);
  	// Function that reads boundary options.
  	bool ReadBoundaryOptions(const char *configFile);
  	// Function that reads temporal options.
  	bool ReadTemporalOptions(const char *configFile);
  	// Function that reads initial condition options.
  	bool ReadICOptions(const char *configFile);
    // Function that reads flow characteristics options.
    bool ReadFlowOptions(const char *configFile);
		// Function that reads solver options.
		bool ReadSolverOptions(const char *configFile);

};
