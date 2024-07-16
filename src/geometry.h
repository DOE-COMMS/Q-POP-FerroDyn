 #pragma once

#include "matrix.h"
#include <iostream>

class geometry_parameters {
private:
	unsigned int* pt_nz_layer_unit = 0;
	long int num_periods; /*number of repeating unit in the superlattice*/
	long int nx_work, ny_work, nz_work;

private:
	void set_geometry();

public:
	bool periodicX, periodicY, periodicZ;
	long int num_layer_unit;
	long int id_lastlayer_unit;
	long int id_firstlayer, id_lastlayer; /*id of the first and the last layer in the superlattice*/

	long int num_layer;
	unsigned int* pt_nz_layer = 0;

	long int nx_system, ny_system, nz_system;

	long int idxi_work = 0, idyi_work = 0, idzi_work = 0;
	long int idxf_work = 0, idyf_work = 0, idzf_work = 0;

	double dx, dy, dz;

	/*PML:
	* Flags to set PML on and off in each direction.
	* Width of the PML is the number of cells in each direction on each side. Usually set to 10.
	*/
	// bool vars set PML at each face (Xs -> YZ plane at X = 0, Xe -> YZ plane at X = Lx)
	bool if_PML, if_PML_Xs, if_PML_Xe, if_PML_Ys, if_PML_Ye, if_PML_Zs, if_PML_Ze;
	// PML parameters
	long int PML_size;
	double sigmaMax; 	// Maximum conductivity in the PML medium
	// Flags for a cell-level check to see if PML is active there
	// matrix3d<bool> isPML;

	long int xS, xE, yS, yE, zS, zE;
	long int nx_phy, ny_phy, nz_phy; // Dimensions of the physical part of the computational mesh

private:
	void set_PML();

public:
	bool if_read_struct;
	static geometry_parameters geo;

public:
	void readgeo();
	void copy_to_device();
	void loggeo();
};
