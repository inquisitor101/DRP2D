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
	// Total number of nodes.
  nNode = geometry_container->GetnNode();

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

	// unsigned short idxElem   = 0;
	// unsigned long  nDOFsGrid = 0;
	// unsigned short nZone     = geometry_container->GetnZone();
  //
  //
	// // Iterate on each zone.
	// for(unsigned short iZone=0; iZone<nZone; iZone++){
  //
	// 	const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);
  //
	// 	unsigned short nPoly     = gridZone->GetnPolySol();
	// 	unsigned short nSubElem  = nPoly*nPoly;
	// 	unsigned long  nElem     = gridZone->GetnElem();
  //
	// 	for(unsigned long iElem=0; iElem<nElem; iElem++){
	// 		const CGeometryElement *surfElem = gridZone->GetGeometryElem(iElem);
	// 		as3data1d<as3double> elemNode = surfElem->GetCoordSolDOFs();
  //
	// 		idxElem = 0;
	// 		for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
	// 			for(unsigned short iNode=0; iNode<4; iNode++){
	// 				for(unsigned short iDim=0; iDim<nDim; iDim++)
	// 					Paraview_File << std::scientific << elemNode[iDim][ConnLocal[iZone][idxElem]] << "\t";
	// 				// accout for z-coordinate
	// 				Paraview_File << std::scientific << "0.0" << "\t";
	// 				idxElem++;
	// 				nDOFsGrid++;
	// 			}
	// 			Paraview_File << "\n";
	// 		}
	// 	}
	// }
  //
  //
	// // Register cell indices.
	// unsigned long nSubElemP1Global 		 = nPoints/N_POINTS_QUADRILATERAL;
	// unsigned long nGlobal_Elem_Storage = nSubElemP1Global*(N_POINTS_QUADRILATERAL+1);
	// Paraview_File << "\nCELLS " << nSubElemP1Global << "\t" << nGlobal_Elem_Storage << "\n";
  //
	// unsigned long idxElemGlobal = 0;
	// for(unsigned short iZone=0; iZone<nZone; iZone++){
	// 	const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);
	// 	unsigned long nElem = gridZone->GetnElem();
  //
	// 	unsigned short nPoly    = gridZone->GetnPolySol();
	// 	unsigned short nSubElem = nPoly*nPoly;
  //
	// 	short idxLocal;
	// 	for(unsigned long iElem=0; iElem<nElem; iElem++){
	// 		idxLocal = 0;
	// 		for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
	// 			Paraview_File << N_POINTS_QUADRILATERAL << "\t";
	// 			for(unsigned short iNode=0; iNode<nNodeP1; iNode++){
	// 				Paraview_File << idxElemGlobal+idxLocal << "\t";
	// 				idxLocal++;
	// 			}
	// 			Paraview_File << "\n";
	// 		}
	// 		idxElemGlobal += idxLocal;
	// 	}
	// }
  //
  //
	// // Cell registration.
	// Paraview_File << "\nCELL_TYPES " << nSubElemP1Global << "\n";
	// for(unsigned short iZone=0; iZone<nZone; iZone++){
	// 	const CGeometryZone *gridZone = geometry_container->GetGeometryZone(iZone);
  //
	// 	unsigned long  nElem 		= gridZone->GetnElem();
	// 	unsigned short nPoly 	  = gridZone->GetnPolySol();
	// 	unsigned short nSubElem = nPoly*nPoly;
  //
	// 	for(unsigned long iElem=0; iElem<nElem; iElem++){
	// 		for(unsigned short iSubElem=0; iSubElem<nSubElem; iSubElem++){
	// 			Paraview_File << QUADRILATERAL << "\t";
	// 			Paraview_File << "\n";
	// 		}
	// 	}
	// }
	// Paraview_File << "\nPOINT_DATA " << nPoints << "\n";
	// std::cout << "Done." << std::endl;
  //
  //
	// // Register density data.
	// std::cout << "  writing Density.........";
	// Paraview_File << "\nSCALARS " << "Density" << " double 1\n";
	// Paraview_File << "LOOKUP_TABLE default\n";
	// WriteScalar(config_container,
	// 						geometry_container,
	// 						solver_container,
	// 						Paraview_File,
	// 						VariableDensity);
	// std::cout << " Done." << std::endl;
  //
	// // Register momentum data.
	// std::cout << "  writing Momentum........";
	// Paraview_File << "\nVECTORS " << "Momentum" << " double\n";
	// WriteVector(config_container,
	// 						geometry_container,
	// 						solver_container,
	// 						Paraview_File,
	// 						VariableMomentum);
	// std::cout << " Done." << std::endl;
  //
	// // Register energy data.
	// std::cout << "  writing Energy..........";
	// Paraview_File << "\nSCALARS " << "Energy" << " double 1\n";
	// Paraview_File << "LOOKUP_TABLE default\n";
	// WriteScalar(config_container,
	// 						geometry_container,
	// 						solver_container,
	// 						Paraview_File,
	// 						VariableEnergy);
	// std::cout << " Done." << std::endl;
  //
	// // Register pressure data.
	// std::cout << "  writing Pressure........";
	// Paraview_File << "\nSCALARS " << "Pressure" << " double 1\n";
	// Paraview_File << "LOOKUP_TABLE default\n";
	// WriteScalar(config_container,
	// 						geometry_container,
	// 						solver_container,
	// 						Paraview_File,
	// 						VariablePressure);
	// std::cout << " Done." << std::endl;
  //
  // // Register temperature data.
	// std::cout << "  writing Temperature.....";
	// Paraview_File << "\nSCALARS " << "Temperature" << " double 1\n";
	// Paraview_File << "LOOKUP_TABLE default\n";
	// WriteScalar(config_container,
	// 						geometry_container,
	// 						solver_container,
	// 						Paraview_File,
	// 						VariableTemperature);
	// std::cout << " Done." << std::endl;
  //
  // // Register Mach number data.
	// std::cout << "  writing Mach............";
	// Paraview_File << "\nSCALARS " << "Mach" << " double 1\n";
	// Paraview_File << "LOOKUP_TABLE default\n";
	// WriteScalar(config_container,
	// 						geometry_container,
	// 						solver_container,
	// 						Paraview_File,
	// 						VariableMachNumber);
	// std::cout << " Done." << std::endl;
  //
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
  //
  //
	// // Close file.
	// Paraview_File.close();
  //
	// // Report: all is complete!
	// std::cout << "Done." << std::endl;
  //
	// // Update file counter.
	// FileNumberVTK++;
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


