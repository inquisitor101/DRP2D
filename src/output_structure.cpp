#include "output_structure.hpp"




COutput::COutput
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor, used to initiate COutput class.
	*/
{
  // Number of nodes in x-dimension.
  nxNode = geometry_container->GetnxNode();
  // Number of nodes in y-dimension.
  nyNode = geometry_container->GetnyNode();
  // Total number of nodes.
  nNode = geometry_container->GetnNode();

  // Number of cells in x-direction.
  nxCell = nxNode-1;
  // Number of cells in y-direction.
  nyCell = nyNode-1;
  // Total number of cells.
  nCell = nxCell*nyCell;

  // Local connectivity mapper according to counter-clockwise convention.
  ConnLocal = { {0, 0}, {1, 0}, {1, 1}, {0, 1} };

	// Extract output VTK filename.
	OutputVTKFilename = config_container->GetOutputVTKFilename();

	// Choose working data to write the output from. Note, for now
	// it is fixed as conservative working-variables.
	VariableDensity     = new CDensityConservative();
	VariableMomentum    = new CMomentumConservative();
	VariableEnergy      = new CEnergyConservative();
	VariablePressure    = new CPressureConservative();
  VariableTemperature = new CTemperatureConservative();
  VariableMachNumber  = new CMachNumberConservative();

	// Reset VTK file number to zero.
	FileNumberVTK = 0;
  // Reset zone data solution file number to zero.
  FileNumberZoneData = 0;
  // Resert processed data file number to zero.
  FileNumberProcessed = 0;
  // Reset data solution file number to zero.
  FileNumberDataSolution = 0;
}


COutput::~COutput
(
 void
)
 /*
	* Destructor for COutput class, frees allocated memory.
	*/
{
	if( VariableDensity     != nullptr ) delete VariableDensity;
	if( VariableMomentum    != nullptr ) delete VariableMomentum;
	if( VariableEnergy      != nullptr ) delete VariableEnergy;
	if( VariablePressure    != nullptr ) delete VariablePressure;
  if( VariableMachNumber  != nullptr ) delete VariableMachNumber;
  if( VariableTemperature != nullptr ) delete VariableTemperature;
}


void COutput::WriteSolutionToFile
(
  CConfig   *config_container,
  CGeometry *geometry_container,
  CSolver   *solver_container,
  as3double  SimTime
)
 /*
  *
  */
{

}


void COutput::WriteDataToFile
(
  CConfig                *config_container,
  CGeometry              *geometry_container,
  const char             *fileinfo,
  as3vector1d<as3double> &time,
  as3vector2d<as3double> &data
)
 /*
  * Function that writes a data profile to a file.
  */
{

}


void COutput::WriteFileVTK
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CSolver   *solver_container
)
 /*
	* Function that writes VTK data file.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution in VTK format... " << std::endl;

	// Create output stream.
	std::ofstream Paraview_File;

	// Open string stream.
	std::stringstream ss;
	ss << OutputVTKFilename << "_" << FileNumberVTK << ".vtk";

	// Open ASCII file.
	Paraview_File.open(ss.str().c_str(), std::ios::out);

	// Check if file can be open.
	if( !Paraview_File.is_open() )
		Terminate("COutput::WriteFileVTK", __FILE__, __LINE__,
							"VTK file directory could not be located!");

	// Write header.
	Paraview_File.precision(6);
	Paraview_File << "# vtk DataFile Version 3.0\n";
	Paraview_File << "vtk output\n";
	Paraview_File << "ASCII\n";
	Paraview_File << "DATASET UNSTRUCTURED_GRID\n";

	// Initialize nPoly=1 element (index) size.
	unsigned short nNodeP1 = N_POINTS_QUADRILATERAL;

	// Extract total number of points assuming nPoly=1 cells.
	unsigned long nPoints = geometry_container->GetnPointCellP1();
	Paraview_File << "POINTS " << nPoints << " double\n";

	// Display progress.
	std::cout << "  writing grid data....... ";

  // Extract coordinates.
  auto* xcoord = geometry_container->GetGridCoordinate(0);
  auto* ycoord = geometry_container->GetGridCoordinate(1);

  // Loop over all cells.
  for(unsigned long j=0; j<nyCell; j++){
    for(unsigned long i=0; i<nxCell; i++){

      // Loop over all nodes in each cell.
      for(unsigned short s=0; s<ConnLocal.size(); s++){

        // Shifted indices to form counter-clockwise convention.
        unsigned long I = i + ConnLocal[s][0];
        unsigned long J = j + ConnLocal[s][1];

        // Global cell node index.
        unsigned long k = J*nxNode+I;
        // Write coordinates.
        Paraview_File << std::scientific << xcoord[k] << "\t"
                                         << ycoord[k] << "\t"
                                         << "0.0"     << "\t";
      }
      Paraview_File << "\n";
    }
  }

  // Register cell indices.
  unsigned long nSubElemP1Global 		 = nPoints/N_POINTS_QUADRILATERAL;
  unsigned long nGlobal_Elem_Storage = nSubElemP1Global*(N_POINTS_QUADRILATERAL+1);
  Paraview_File << "\nCELLS " << nSubElemP1Global << "\t" << nGlobal_Elem_Storage << "\n";

  // Loop over all cells.
  unsigned long idxElemGlobal = 0;
  for(unsigned long i=0; i<nCell; i++){
    // Register number of points needed in each cell.
    Paraview_File << N_POINTS_QUADRILATERAL << "\t";
    // Loop over all nodes in each cell.
    for(unsigned short j=0; j<ConnLocal.size(); j++){
        // Register local cell-node index convention.
        Paraview_File << idxElemGlobal + j << "\t";
    }
		Paraview_File << "\n";
    // Update global index.
    idxElemGlobal += ConnLocal.size();
  }

  // Cell registration.
	Paraview_File << "\nCELL_TYPES " << nSubElemP1Global << "\n";
  // Loop over all cells.
  for(unsigned long i=0; i<nCell; i++){
    // Register number of points needed in each cell.
    Paraview_File << QUADRILATERAL << "\n";
  }
  Paraview_File << "\nPOINT_DATA " << nPoints << "\n";
  std::cout << "Done." << std::endl;


	// Register density data.
	std::cout << "  writing Density.........";
	Paraview_File << "\nSCALARS " << "Density" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableDensity);
	std::cout << " Done." << std::endl;


	// Register momentum data.
	std::cout << "  writing Momentum........";
	Paraview_File << "\nVECTORS " << "Momentum" << " double\n";
	WriteVector(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableMomentum);
	std::cout << " Done." << std::endl;

	// Register energy data.
	std::cout << "  writing Energy..........";
	Paraview_File << "\nSCALARS " << "Energy" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableEnergy);
	std::cout << " Done." << std::endl;

	// Register pressure data.
	std::cout << "  writing Pressure........";
	Paraview_File << "\nSCALARS " << "Pressure" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariablePressure);
	std::cout << " Done." << std::endl;

  // Register temperature data.
	std::cout << "  writing Temperature.....";
	Paraview_File << "\nSCALARS " << "Temperature" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableTemperature);
	std::cout << " Done." << std::endl;

  // Register Mach number data.
	std::cout << "  writing Mach............";
	Paraview_File << "\nSCALARS " << "Mach" << " double 1\n";
	Paraview_File << "LOOKUP_TABLE default\n";
	WriteScalar(config_container,
							geometry_container,
							solver_container,
							Paraview_File,
							VariableMachNumber);
	std::cout << " Done." << std::endl;
  
  // // If a PML zone exists, write its auxiliary variables.
  // if( config_container->GetUsePML() ){
  //
  //   // If user specifies auxiliary data included, write them.
  //   if( config_container->GetWriteAuxiliaryDataPML() ){
  //
  //     // Register Q1 of density data.
  //     std::cout << "  writing Q1[Density].....";
  //     Paraview_File << "\nSCALARS " << "Q1[Density]" << " double 1\n";
  //     Paraview_File << "LOOKUP_TABLE default\n";
  //     WriteScalarAuxPML(config_container,
  //                       geometry_container,
  //                       solver_container,
  //                       Paraview_File,
  //                       CONTQ1_VAR);
  //     std::cout << " Done." << std::endl;
  //
  //     // Register Q1 of momentum data.
  //     std::cout << "  writing Q1[Momentum]....";
  //     Paraview_File << "\nVECTORS " << "Q1[Momentum]" << " double\n";
  //     WriteVectorAuxPML(config_container,
  //                       geometry_container,
  //                       solver_container,
  //                       Paraview_File,
  //                       XMOMQ1_VAR);
  //     std::cout << " Done." << std::endl;
  //
  //     // Register Q1 of energy data.
  //     std::cout << "  writing Q1[Energy]......";
  //     Paraview_File << "\nSCALARS " << "Q1[Energy]" << " double 1\n";
  //     Paraview_File << "LOOKUP_TABLE default\n";
  //     WriteScalarAuxPML(config_container,
  //                       geometry_container,
  //                       solver_container,
  //                       Paraview_File,
  //                       ENERQ1_VAR);
  //     std::cout << " Done." << std::endl;
  //   }
  // }


	// Close file.
	Paraview_File.close();

	// Report: all is complete!
	std::cout << "Done." << std::endl;

	// Update file counter.
	FileNumberVTK++;
}


void COutput::WriteScalar
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		   *solver_container,
 std::ofstream   &Paraview_File,
 CScalarVariable *Variable
)
 /*
	* Function that writes a scalar parameter.
	*/
{
  // Extract solution data.
  const auto& variables = solver_container->GetDataSolution();

  // Loop over all cells.
  for(unsigned long j=0; j<nyCell; j++){
    for(unsigned long i=0; i<nxCell; i++){

      // Loop over all nodes in each cell.
      for(unsigned short s=0; s<ConnLocal.size(); s++){

        // Shifted indices to form counter-clockwise convention.
        unsigned long I = i + ConnLocal[s][0];
        unsigned long J = j + ConnLocal[s][1];

        // Global cell node index.
        unsigned long k = J*nxNode+I;

        // Compute scalar.
        auto scalar = Variable->GetValue(variables, k);
        // Write coordinates.
        Paraview_File << std::scientific << scalar << "\t";
      }
      Paraview_File << "\n";
    }
  }
}


void COutput::WriteVector
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		   *solver_container,
 std::ofstream   &Paraview_File,
 CVectorVariable *Variable
)
 /*
	* Function that writes a vector parameter.
	*/
{
  // Extract solution data.
  const auto& variables = solver_container->GetDataSolution();

  // Loop over all cells.
  for(unsigned long j=0; j<nyCell; j++){
    for(unsigned long i=0; i<nxCell; i++){

      // Loop over all nodes in each cell.
      for(unsigned short s=0; s<ConnLocal.size(); s++){

        // Shifted indices to form counter-clockwise convention.
        unsigned long I = i + ConnLocal[s][0];
        unsigned long J = j + ConnLocal[s][1];

        // Global cell node index.
        unsigned long k = J*nxNode+I;

        // Compute scalar.
        auto vector = Variable->GetValue(variables, k);
        // Write coordinates.
        Paraview_File << std::scientific << vector[0] << "\t"
                                         << vector[1] << "\t"
                                         << "0.0"     << "\t";
      }
      Paraview_File << "\n";
    }
  }
}


void COutput::WriteScalarAuxPML
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		   *solver_container,
 std::ofstream   &Paraview_File,
 unsigned short   iVar
)
 /*
	* Function that writes a PML auxiliary scalar parameter.
	*/
{

}


void COutput::WriteVectorAuxPML
(
 CConfig   		   *config_container,
 CGeometry 		   *geometry_container,
 CSolver  		   *solver_container,
 std::ofstream   &Paraview_File,
 unsigned short   iVar
)
 /*
	* Function that writes a PML auxiliary vector parameter.
	*/
{

}


