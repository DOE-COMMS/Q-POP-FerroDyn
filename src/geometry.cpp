#include "geometry.h"

void geometry_parameters::set_geometry() {
	long int num_unitlayers;

	num_unitlayers = id_lastlayer_unit - id_firstlayer + 1;
	id_lastlayer = id_firstlayer - 1 + num_periods * num_unitlayers;
	num_layer = num_layer_unit + (num_periods - 1) * num_unitlayers;

	pt_nz_layer = new unsigned int[num_layer];
	long int j = id_firstlayer;
	for (long int i = 0; i < num_layer; i++) {
		if (i < id_firstlayer - 1) {
			pt_nz_layer[i] = pt_nz_layer_unit[i];
		}
		else if (i > id_lastlayer - 1) {
			pt_nz_layer[i] = pt_nz_layer_unit[i - num_layer + num_layer_unit];
		}
		else {
			pt_nz_layer[i] = pt_nz_layer_unit[j - 1];
			j = j + 1;
			if (j > id_lastlayer_unit) {
				j = id_firstlayer;
			}
		}
	}

	nz_work = 0;
	for (long int i = 0; i < num_layer; i++) {
		nz_work = nz_work + pt_nz_layer[i];
	}

	if (nz_work > nz_system) {
		std::cout << "Layer thickness sums to " << nz_work << ", which exceeds length of the system nz=" << nz_system << std::endl;
		exit(1);
	}

	set_PML();

	idxi_work = static_cast<long int>((nx_system - nx_work) / 2) + 1;
	idxf_work = idxi_work + nx_work - 1;
	idyi_work = static_cast<long int>((ny_system - ny_work) / 2) + 1;
	idyf_work = idyi_work + ny_work - 1;
	//idzi_work = static_cast<long int>((nz_system - nz_work) / 2) + 1;
	//idzf_work = idzi_work + nz_work - 1;
	idzi_work = (if_PML_Zs? PML_size : 0) + 1;
	idzf_work = idzi_work + nz_work - 1;

	idxi_work = idxi_work - 1; idyi_work = idyi_work - 1; idzi_work = idzi_work - 1;
	idxf_work = idxf_work - 1; idyf_work = idyf_work - 1; idzf_work = idzf_work - 1;
}

void geometry_parameters::set_PML()
{
	if (if_PML == false)
	{
		PML_size = 0;
		PML_materialType = 0;

		if_PML_Xs = false;
		if_PML_Xe = false;
		if_PML_Ys = false;
		if_PML_Ye = false;
		if_PML_Zs = false;
		if_PML_Ze = false;

		xS = 0;
		xE = nx_system;

		yS = 0;
		yE = ny_system;

		zS = 0;
		zE = nz_system;
	}
	else
	{
		if (periodicX && (if_PML_Xs || if_PML_Xe))
		{
			std::cout << "Cannot have PML and periodicity both in X direction. Defaulting to periodic." << std::endl;
			if_PML_Xs = false;
			if_PML_Xe = false;
		}

		if (periodicY && (if_PML_Ys || if_PML_Ye))
		{
			std::cout << "Cannot have PML and periodicity both in Y direction. Defaulting to periodic." << std::endl;
			if_PML_Ys = false;
			if_PML_Ye = false;
		}

		if (periodicZ && (if_PML_Zs || if_PML_Ze))
		{
			std::cout << "Cannot have PML and periodicity both in Z direction. Defaulting to periodic." << std::endl;
			if_PML_Zs = false;
			if_PML_Ze = false;
		}

		if (if_PML_Xs == false && if_PML_Xe == false \
			&& if_PML_Ys == false && if_PML_Ye == false \
			&& if_PML_Zs == false && if_PML_Ze == false)
		{
			if_PML = false;
		}

		if (if_PML_Xs || if_PML_Xe)
		{
			if_PML_Xs = true;
			if_PML_Xe = true;
		}

		if (if_PML_Ys || if_PML_Ye)
		{
			if_PML_Ys = true;
			if_PML_Ye = true;
		}

		if (if_PML_Zs || if_PML_Ze)
		{
			if_PML_Zs = true;
			if_PML_Ze = true;
		}

		xS = (if_PML_Xs ? PML_size : 0);
		xE = (if_PML_Xs ? nx_system + PML_size : nx_system);

		yS = (if_PML_Ys ? PML_size : 0);
		yE = (if_PML_Ys ? ny_system + PML_size : ny_system);

		zS = (if_PML_Zs ? PML_size : 0);
		zE = (if_PML_Zs ? nz_system + PML_size : nz_system);

		if (if_PML_Xs)
			nx_system += PML_size;
		if (if_PML_Xe)
			nx_system += PML_size;

		if (if_PML_Ys)
			ny_system += PML_size;
		if (if_PML_Ye)
			ny_system += PML_size;

		if (if_PML_Zs)
			nz_system += PML_size;
		if (if_PML_Ze)
			nz_system += PML_size;
	}

	nx_phy = nx_system - (if_PML_Xs ? PML_size : 0) - (if_PML_Xe ? PML_size : 0);
	ny_phy = ny_system - (if_PML_Ys ? PML_size : 0) - (if_PML_Ye ? PML_size : 0);
	nz_phy = nz_system - (if_PML_Zs ? PML_size : 0) - (if_PML_Ze ? PML_size : 0);
}

void geometry_parameters::copy_to_device() {
#pragma acc enter data copyin(this)
}

geometry_parameters geometry_parameters::geo;
