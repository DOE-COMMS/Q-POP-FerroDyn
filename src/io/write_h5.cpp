#include "write_h5.h"


/*
 *  1. matrix3d data
 *  2. Time step to write.
 *  3, 4, 5. Mesh size
 *  6, 7, 8. Cell size
 *  9. File root (no timestep or filetype extension)
 *  10. Variable name
 *  13. appendFlag = 0 -> write a new file, appendFlag != 0 -> add field to existing file
 */
bool saveVTK_Scalar(matrix3d<double> data,
    int t,
	long int xS, long int yS, long int zS,
	long int xE, long int yE, long int zE,
	double dx, double dy, double dz,
    const std::string& fileRoot, const std::string& varName,
    int appendFlag)
{
	long int nx = xE - xS;
	long int ny = yE - yS;
	long int nz = zE - zS;

	// Validate input parameters
    if (nx <= 0 || ny <= 0 || nz <= 0 || dx <= 0 || dy <= 0 || dz <= 0 || fileRoot.empty() || varName.empty())
    {
        std::cerr << "Invalid input parameters." << std::endl;
        return false;
    }

    std::string fileName = fileRoot + '_' + std::to_string(t) + ".vtk";

    std::ofstream vtkFile;
    if (appendFlag == 0)
    {
        vtkFile.open(fileName.c_str(), std::ios::out);
    }
    else
    {
        vtkFile.open(fileName.c_str(), std::ios::app);
    }

    if (!vtkFile.is_open())
    {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return false;
    }
    else
    {
        std::cout << "Writing " << fileName << std::endl;
    }

    // VTK header
    if (appendFlag == 0)
    {
        vtkFile << "# vtk DataFile Version 3.0" << std::endl;
        vtkFile << "VTK file for time step " << t << std::endl;
        vtkFile << "ASCII" << std::endl;
        vtkFile << "DATASET STRUCTURED_POINTS" << std::endl;
        vtkFile << "DIMENSIONS " << ny << " " << nx << " " << nz << std::endl;
        vtkFile << "ORIGIN 0 0 0" << std::endl;
        vtkFile << "SPACING " << dy << " " << dx << " " << dz << std::endl;
        vtkFile << "POINT_DATA " << nx * ny * nz << std::endl;
    }

    vtkFile << std::endl << "SCALARS " << varName << " double 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				vtkFile << data(xS + i, yS + j, zS + k) << '\n';
			}
		}
	}

    vtkFile.close();
    return true;
}

bool saveVTK_Vector(matrix3d<double> data1, matrix3d<double> data2, matrix3d<double> data3,
    int t,
	long int xS, long int yS, long int zS,
	long int xE, long int yE, long int zE,
	double dx, double dy, double dz,
    const std::string& fileRoot, const std::string& varName,
    int appendFlag)
{
	long int nx = xE - xS;
	long int ny = yE - yS;
	long int nz = zE - zS;

	// Validate input parameters
    if (nx <= 0 || ny <= 0 || nz <= 0 || dx <= 0 || dy <= 0 || dz <= 0 || fileRoot.empty() || varName.empty())
    {
        std::cerr << "Invalid input parameters." << std::endl;
        return false;
    }
    std::string fileName = fileRoot + '_' + std::to_string(t) + ".h5";

	using namespace HighFive;

    File file(fileName, File::ReadWrite | File::Create | File::Truncate);

	std::vector<double> data(nz);

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				data[k] = data1(xS + i, yS + j, zS + k);
			}
		}
	}
	
	DataSet dataset = file.createDataSet<double>(varName, DataSpace::From(data));
	dataset.write(data);

    return true;
}