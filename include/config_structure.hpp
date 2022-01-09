#pragma once

#include "option_structure.hpp"
#include <map>



class CConfig {

	public:
		// Constructor.
		CConfig(const char *configFile);

		// Destructor.
		~CConfig(void);

    // Getter: returns nZone.
    unsigned short GetnZone(void)                              const {return nZone;}
    // Getter: returns nxNode[iZone].
    unsigned long GetnxNode(unsigned short iZone)              const {return nxNode[iZone];}
    // Getter: returns nyNode[iZone].
    unsigned long GetnyNode(unsigned short iZone)              const {return nyNode[iZone];}
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
    // Getter: returns TypeSolver[iZone].
    unsigned short GetTypeSolver(unsigned short iZone)         const {return TypeSolver[iZone];}
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
    // Getter: returns UniformGridResolution.
    bool GetUniformGridResolution(void)                        const {return UniformGridResolution;}
    // Getter: returns ZoneMarker, per input zone.
    unsigned short GetTypeZone(unsigned short iZone)           const {return TypeZone[iZone];}
    // Getter: returns dsNodalRatioZone.
    const as3vector1d<as3double> GetdsNodalRatioZone(void)     const {return dsNodalRatioZone;}

	private:
    // Default data parameters.
		struct {
      // Number of zones.
      unsigned short nZone = 1;
			// Zone marker.
			as3vector1d<std::string> NameZoneMarker = { "ZONE_MAIN" };
			// Markers.
			as3vector1d<std::string> NameMarker = {};
			// Restart solution or not.
			std::string NameRestartSolution = "false";
      // Ratio of nodal sizes in (WEST, EAST, SOUTH, NORTH).
      as3vector1d<as3double> dsNodalRatioZone = { 1.0, 1.0, 1.0, 1.0 };
      // Zone conformity w.r.t. ZONE_MAIN.
      std::string NameZoneConformity = "true";
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
      // Uniform equidistant grid resolution or not.
      std::string NameUniformGridResolution = "true";

		} DefaultParam;

    // Number of zones.
    unsigned short nZone;
    // Polynomial order of solution per zone.
    as3vector1d<unsigned short> nPolySolZone;
    // Multizone strategy used in the simulation.
    unsigned short MultizoneStrategy;

    // Zone markers, string name.
    as3vector1d<std::string> NameZoneMarker;
    // Zone type, enum values.
    as3vector1d<unsigned short> TypeZone;

		// Number of nodes in x-direction, per zone.
		as3vector1d<unsigned long> nxNode;
		// Number of nodes in y-direction, per zone.
		as3vector1d<unsigned long> nyNode;

    // Nodal ratio per extension of zone: (WEST, EAST, SOUTH, NORTH).
    as3vector1d<as3double> dsNodalRatioZone;
    // Use zone conformity between adjacent nodes on zone.
    bool ZoneConformity;
    // Uniform equidistant grid resolution or not.
    bool UniformGridResolution;

    // Type of stencils used in eaech direction, string name.
    as3vector1d<std::string> NameTypeStencil;
    // Type of stencil used in x-direction.
    as3vector1d<unsigned short> TypeStencil;

    // Domain bounds (physical domain).
    // Convention: (WEST, EAST, SOUTH, NORTH).
    as3vector1d<as3double> DomainBound;
    // Type of solver per zone.
    as3vector1d<unsigned short> TypeSolver;

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


    // Function that maps TypeZone from string to enum.
    void MapTypeZone(void);
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


    // Function that determines the multizone strategy adopted.
    void DetermineMultizoneStrategy(void);
    // Function that checks the nodal ratios specified according
    // to specified multizone strategy.
    void CheckdsNodalRatio(void);
    // Function that processes the zone conformity in the specified grids.
    void ProcessZoneConformity(void);


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
