#include "global.h"

void global_parameters::set_global() {
	long int id_zi, id_zf;
	long int ti, tj, tk;
	unsigned int mat_type, mat_typei, mat_typef;
	unsigned long long ct;
	material* mat, * mati, * matf;
	bool if_spinsource_i, if_spinsource_f;
	double norm_dm;

	if (false == if_EMdynamic) pt_geo->if_PML = false;

	pt_geo->set_PML();
	
	nx = pt_geo->nx_system;
	ny = pt_geo->ny_system;
	nz = pt_geo->nz_system;

	nx_phy = pt_geo->nx_phy;
	ny_phy = pt_geo->ny_phy;
	nz_phy = pt_geo->nz_phy;

	n = nx * ny * nz;

	dx = pt_geo->dx;
	dy = pt_geo->dy;
	dz = pt_geo->dz;

	if (pt_geo->if_PML) {
		kappaMax = pt_geo->kappaMax;
		PML_m = pt_geo->PML_m;

		if (pt_geo->PML_materialType > 0 && pt_geo->PML_materialType < num_materials)
		{
			mat = &(material_parameters[pt_geo->PML_materialType - 1]);
			PML_er11 = mat->r_permittivity11;
			PML_er22 = mat->r_permittivity22;
			PML_er33 = mat->r_permittivity33;
		}
		else
		{
			PML_er11 = 1.;
			PML_er22 = 1.;
			PML_er33 = 1.;
		}

		eta0 = sqrt(mu0 / e0 * ((PML_er11 + PML_er22 + PML_er33) / 3.0));
		maxReflErr = exp(-16);

		if (pt_geo->PML_size != 0)
			sigmaMax = -(PML_m + 1.0) * log(maxReflErr) / (2.0 * pt_geo->PML_size * eta0 * (dx + dy + dz) / 3.0);
		else
			sigmaMax = 0.0;

		if (pt_geo->if_PML_Xs)
		{
			Jfin_xi += pt_geo->PML_size;
			Jfin_xf += pt_geo->PML_size;
		}

		if (pt_geo->if_PML_Ys)
		{
			Jfin_yi += pt_geo->PML_size;
			Jfin_yf += pt_geo->PML_size;
		}

		if (pt_geo->if_PML_Zs)
		{
			Jfin_zi += pt_geo->PML_size;
			Jfin_zf += pt_geo->PML_size;
		}
	}

	Hext[0] = 0.; Hext[1] = 0.; Hext[2] = 0.;
	Eext[0] = 0.; Eext[1] = 0.; Eext[2] = 0.;

	//---------Designate Material Type for each layer-----------//
	pt_material_layer = new unsigned int[pt_geo->num_layer];
	long int j = pt_geo->id_firstlayer;
	for (long int i = 0; i < pt_geo->num_layer; i++) {
		if (i < pt_geo->id_firstlayer - 1) {
			pt_material_layer[i] = pt_material_layer_unit[i];
		}
		else if (i > pt_geo->id_lastlayer - 1) {
			pt_material_layer[i] = pt_material_layer_unit[i - pt_geo->num_layer + pt_geo->num_layer_unit];
		}
		else {
			pt_material_layer[i] = pt_material_layer_unit[j - 1];
			j = j + 1;
			if (j > pt_geo->id_lastlayer_unit) {
				j = pt_geo->id_firstlayer;
			}
		}
	}
	//---------END Designate Material Type for each layer-----------//
	// 
	//---------Designate Material Type for each cell-----------//
	material_cell.initialize(nx, ny, nz);

	if (pt_geo->if_read_struct == false) {
		id_zf = pt_geo->idzi_work - 1;
		for (long int id = 0; id < pt_geo->num_layer; id++) {
			id_zi = id_zf + 1; id_zf = id_zi + pt_geo->pt_nz_layer[id] - 1;
			for (long int i = pt_geo->idxi_work; i < pt_geo->idxf_work + 1; i++) {
				for (long int j = pt_geo->idyi_work; j < pt_geo->idyf_work + 1; j++) {
					for (long int k = id_zi; k < id_zf + 1; k++) {
						material_cell(i, j, k) = pt_material_layer[id];
					}
				}
			}
		}
	}
	else {
		read_struct();
	}
	//---------END Designate Material Type for each cell-----------//

	//-------Determine if there is FM(E) and AFM in system------//
	if_FM_all = false; if_FE_all = false; if_AFM_all = false;
	for (long int i = 0; i < num_materials; i++) {
		if (material_parameters[i].if_FM == true) {
			if_FM_all = true;
		}
		if (material_parameters[i].if_FE == true) {
			if_FE_all = true;
		}
		if (material_parameters[i].if_AFM == true) {
			if_AFM_all = true;
		}
	}

	NFE = 0; NFM = 0; NAFM = 0;
	for (long id = 0; id < n; id++) {
		mat_type = material_cell(id);
		if (mat_type != 0) {
			mat = &(material_parameters[mat_type - 1]);
			if (mat->if_FM == true) {
				NFM = NFM + 1;
			}
			if (mat->if_FE == true) {
				NFE = NFE + 1;
			}
			if (mat->if_AFM == true) {
				NAFM = NAFM + 1;
			}
		}
	}

	if (if_output_only_magcell == true) {
		magcell_index = new unsigned long long[NFM];
		ct = 0;
		for (long id = 0; id < n; id++) {
			mat_type = material_cell(id);
			if (mat_type != 0) {
				mat = &(material_parameters[mat_type - 1]);
				if (mat->if_FM == true) {
					ct = ct + 1;
					magcell_index[ct - 1] = id;
				}
			}
		}
	}

	//----------Determine if all surfaces have periodic conditions---//
	if_periodic_allsurface = false;
	if (pt_geo->periodicX == true && pt_geo->periodicY == true && pt_geo->periodicZ == true) {
		if_periodic_allsurface = true;
	}

	/* Prioritize PML over ABC */
	if (pt_geo->if_PML == true && if_EMdynamic == true)
	{
		if_1D_ABC = false;
		if_1D_ABC_ELAST_onlytop = false;
		if_1D_ABC_EM_onlytop = false;
	}

	//---------Initialization for surface free charge density----------//
	if(pt_geo->if_PML_Zs && pt_geo->if_PML) {
		free_charge_nzi += pt_geo->PML_size;
		free_charge_nzf += pt_geo->PML_size;
	}
	free_charge_surfaces = new double[nz + 1];
	free_charge_cells = new double[nz];
	for (long i = 0; i < nz + 1; i++) {
		if (i >= free_charge_nzi && i < free_charge_nzf) {
			free_charge_surfaces[i] = free_charge;
		}
		else {
			free_charge_surfaces[i] = 0.;
		}
	}

	for (long i = 0; i < nz; i++) {
		free_charge_cells[i] = 0.;
		for (long j = 0; j < i+1; j++) {
			free_charge_cells[i] = free_charge_cells[i] + free_charge_surfaces[j];
		}
		for (long j = i+1; j < nz + 1; j++) {
			free_charge_cells[i] = free_charge_cells[i] - free_charge_surfaces[j];
		}
	}

	//----------Initialization for spin pumping and inverse spin Hall effect(ISHE) induced current
	if (if_spin_pumping == true) {
		if_J_ISHE = new bool[nz];
		source_z = new long[nz];
		thickness_pumping_layer = new double[nz];

		for (unsigned int i = 0; i < nz; i++) {
			if_J_ISHE[i] = false;
			source_z[i] = -1;
			thickness_pumping_layer[i] = 0.;
		}

		id_zf = pt_geo->idzi_work - 1;
		for (unsigned int id = 0; id < pt_geo->num_layer; id++) {
			id_zi = id_zf + 1; id_zf = id_zi + pt_geo->pt_nz_layer[id] - 1;
			mat_type = material_cell(pt_geo->if_PML_Xs? pt_geo->PML_size:0, pt_geo->if_PML_Ys? pt_geo->PML_size:0, id_zi);
			mat = &(material_parameters[(mat_type)-1]);
			if (mat->if_spin_pump == true) {
				if (id_zi != (pt_geo->if_PML_Zs? pt_geo->PML_size:0)) {
					mat_typei = material_cell(pt_geo->if_PML_Xs? pt_geo->PML_size:0, pt_geo->if_PML_Ys? pt_geo->PML_size:0, id_zi - 1);
					if (mat_typei != 0) {
						mati = &(material_parameters[(mat_typei)-1]);
						if (mati->if_FM == true || mati->if_AFM == true) {
							if_spinsource_i = true;
						}
						else {
							if_spinsource_i = false;
						}
					}
					else {
						if_spinsource_i = false;
					}
				}
				else {
					if_spinsource_i = false;
				}

				if (id_zf != nz - (pt_geo->if_PML_Ze? pt_geo->PML_size+1:1)) {
					mat_typef = material_cell(pt_geo->if_PML_Xs? pt_geo->PML_size:0, pt_geo->if_PML_Ys? pt_geo->PML_size:0, id_zf + 1);
					if (mat_typef != 0) {
						matf = &(material_parameters[(mat_typef)-1]);
						if (matf->if_FM == true || matf->if_AFM == true) {
							if_spinsource_f = true;
						}
						else {
							if_spinsource_f = false;
						}
					}
					else {
						if_spinsource_f = false;
					}
				}
				else {
					if_spinsource_f = false;
				}

				if (if_spinsource_i == true && if_spinsource_f == false) {
					for (long k = id_zi; k < id_zf + 1; k++) {
						if_J_ISHE[k] = true;
						source_z[k] = id_zi - 1;
						thickness_pumping_layer[k] = pt_geo->pt_nz_layer[id] * dz;
					}
				}
				else if (if_spinsource_i == false && if_spinsource_f == true) {
					for (long k = id_zi; k < id_zf + 1; k++) {
						if_J_ISHE[k] = true;
						source_z[k] = id_zf + 1;
						thickness_pumping_layer[k] = pt_geo->pt_nz_layer[id] * dz;
					}
				}

			}
		}
	}

	//-------Non-uniform external magnetic field-----------//
	Hextx_nonunif.initialize(nx, ny, nz);
	Hexty_nonunif.initialize(nx, ny, nz);
	Hextz_nonunif.initialize(nx, ny, nz);

	if (if_read_Hext_nonunif == true) {
		read_Hext_nonunif();
	}

	//-------Prescribed external electric field-----------//
	if (if_prescribe_Eext == true) {
		prescribe_Eext_time = new double[num_prescribe_Eext];
		prescribe_Ex = new double[num_prescribe_Eext];
		prescribe_Ey = new double[num_prescribe_Eext];
		prescribe_Ez = new double[num_prescribe_Eext];

		read_prescribe_Eext();
	}
	//-------Prescribed external magnetic field-----------//
	if (if_prescribe_Hext == true) {
		prescribe_Hext_time = new double[num_prescribe_Hext];
		prescribe_Hx = new double[num_prescribe_Hext];
		prescribe_Hy = new double[num_prescribe_Hext];
		prescribe_Hz = new double[num_prescribe_Hext];

		read_prescribe_Hext();
	}

	//-------Normalize demag factors------------------------//
	norm_dm = demag_fac_x + demag_fac_y + demag_fac_z;
	demag_fac_x = demag_fac_x / norm_dm * -1.;
	demag_fac_y = demag_fac_y / norm_dm * -1.;
	demag_fac_z = demag_fac_z / norm_dm * -1.;

	//-------Initialize variables in frequency domain-------//
	dkx = 2. * PI / (static_cast<double>(nx) * dx);
	dky = 2. * PI / (static_cast<double>(ny) * dy);
	dkz = 2. * PI / (static_cast<double>(nz) * dz);

	kx = new double[nx];
	ky = new double[ny];
	kz = new double[nz / 2 + 1];

	for (long i = 1; i < nx + 1; i++) {
		ti = (i)-1;
		if (ti > static_cast<long>(nx) / 2) {
			ti = ti - static_cast<long>(nx);
		}
		kx[i - 1] = static_cast<double>(ti) * dkx;
	}

	for (long j = 1; j < ny + 1; j++) {
		tj = (j)-1;
		if (tj > static_cast<long>(ny) / 2) {
			tj = tj - static_cast<long>(ny);
		}
		ky[j - 1] = static_cast<double>(tj) * dky;
	}

	for (long k = 1; k < nz / 2 + 2; k++) {
		tk = (k)-1;
		if (tk > static_cast<long>(nz) / 2) {
			tk = tk - static_cast<long>(nz);
		}
		kz[k - 1] = static_cast<double>(tk) * dkz;
	}

	k_norm.initialize(nx, ny, nz / 2 + 1);

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz / 2 + 1; k++) {
				k_norm(i, j, k) = pow(kz[k], 2.) + pow(ky[j], 2.) + pow(kx[i], 2.);
			}
		}
	}
}

void global_parameters::check_material(material& mat) {
	if (mat.if_FM == false) {
		mat.Ms = 0.; mat.gyro = 0.; mat.FM_damping = 0.;
		mat.anisotropy_type = 0;
		mat.K1 = 0.; mat.K2 = 0.;
		mat.Aex = 0.;
		mat.iDMI = 0.;
		mat.lamda100 = 0.; mat.lamda111 = 0.;
		mat.B1 = 0.; mat.B2 = 0.;
	}
	if (mat.if_AFM == false) {
		mat.Ms_AFM1 = 0.; mat.Ms_AFM2 = 0.;
		mat.gyro_AFM1 = 0.; mat.gyro_AFM2 = 0.;
		mat.damping_AFM1 = 0.; mat.damping_AFM2 = 0.;
		mat.K1_AFM1 = 0.; mat.K1_AFM2 = 0.;
		mat.Aex_AFM1 = 0.; mat.Aex_AFM2 = 0.;
		mat.J_AFM = 0.;
		mat.lamda100_AFM1 = 0.; mat.lamda100_AFM2 = 0.;
		mat.lamda111_AFM1 = 0.; mat.lamda111_AFM2 = 0.;
		mat.B1_AFM1 = 0.; mat.B1_AFM2 = 0.;
		mat.B2_AFM1 = 0.; mat.B2_AFM2 = 0.;
	}
	if (mat.if_FE == false) {
		mat.FE_mass = 0.; mat.FE_damping = 0.;
		mat.a1 = 0.;
		mat.a11 = 0.;  mat.a12 = 0.;
		mat.a111 = 0.;  mat.a112 = 0.;  mat.a123 = 0.;
		mat.a1111 = 0.;  mat.a1112 = 0.;  mat.a1122 = 0.;  mat.a1123 = 0.;
		mat.G11 = 0.;
		mat.Q11 = 0.; mat.Q12 = 0.;  mat.Q44 = 0.;

		mat.tv1 = 0.; mat.tv2 = 0.; mat.tv3 = 0.; mat.tv4 = 0.;
		mat.tv5 = 0.; mat.tv6 = 0.; mat.tv7 = 0.; mat.tv8 = 0.;

		mat.f11 = 0.; mat.f12 = 0.;  mat.f44 = 0.;

		mat.F11 = 0.; mat.F12 = 0.; mat.F44 = 0.;
	}

	mat.comp_n1 = mat.comp_n1 / (mat.comp_n1 + mat.comp_n2 + mat.comp_n3);
	mat.comp_n2 = mat.comp_n2 / (mat.comp_n1 + mat.comp_n2 + mat.comp_n3);
	mat.comp_n3 = mat.comp_n3 / (mat.comp_n1 + mat.comp_n2 + mat.comp_n3);
}

void global_parameters::get_cijkl(material& mat) {
	double Ccrt[6][6] = { 0. };
	double Cijkl[3][3][3][3] = { 0. }, Cwxyz[3][3][3][3] = { 0. };
	unsigned int i, j, k, l, w, x, y, z;

	Ccrt[0][0] = mat.c11; Ccrt[1][1] = mat.c11; Ccrt[2][2] = mat.c11;
	Ccrt[0][1] = mat.c12; Ccrt[0][2] = mat.c12; Ccrt[1][2] = mat.c12;
	Ccrt[1][0] = mat.c12; Ccrt[2][0] = mat.c12; Ccrt[2][1] = mat.c12;
	Ccrt[3][3] = mat.c44; Ccrt[4][4] = mat.c44; Ccrt[5][5] = mat.c44;

	for (i = 1; i < 7; i++) {
		for (j = 1; j < 7; j++) {
			switch (i) {
			case 1:	w = 1; x = 1; break;
			case 2:	w = 2; x = 2; break;
			case 3:	w = 3; x = 3; break;
			case 4:	w = 2; x = 3; break;
			case 5:	w = 1; x = 3; break;
			case 6:	w = 1; x = 2; break;
			}
			switch (j) {
			case 1:	y = 1; z = 1; break;
			case 2:	y = 2; z = 2; break;
			case 3: y = 3; z = 3; break;
			case 4:	y = 2; z = 3; break;
			case 5:	y = 1; z = 3; break;
			case 6:	y = 1; z = 2; break;
			}
			Cijkl[w - 1][x - 1][y - 1][z - 1] = Ccrt[i - 1][j - 1];
			Cijkl[x - 1][w - 1][y - 1][z - 1] = Ccrt[i - 1][j - 1];
			Cijkl[w - 1][x - 1][z - 1][y - 1] = Ccrt[i - 1][j - 1];
			Cijkl[x - 1][w - 1][z - 1][y - 1] = Ccrt[i - 1][j - 1];
		}
	}
	for (w = 0; w < 3; w++) {
		for (x = 0; x < 3; x++) {
			for (y = 0; y < 3; y++) {
				for (z = 0; z < 3; z++) {
					for (i = 0; i < 3; i++) {
						for (j = 0; j < 3; j++) {
							for (k = 0; k < 3; k++) {
								for (l = 0; l < 3; l++) {
									Cwxyz[w][x][y][z] = Cwxyz[w][x][y][z] + \
										Tc2g[w][i] * Tc2g[x][j] * Tc2g[y][k] * Tc2g[z][l] \
										* Cijkl[i][j][k][l];
								}
							}
						}
					}
				}
			}
		}
	}

	for (i = 1; i < 7; i++) {
		for (j = 1; j < 7; j++) {
			switch (i) {
			case(1):w = 1; x = 1; break;
			case(2):w = 2; x = 2; break;
			case(3):w = 3; x = 3; break;
			case(4):w = 2; x = 3; break;
			case(5):w = 1; x = 3; break;
			case(6):w = 1; x = 2; break;
			}

			switch (j) {
			case(1):y = 1; z = 1; break;
			case(2):y = 2; z = 2; break;
			case(3):y = 3; z = 3; break;
			case(4):y = 2; z = 3; break;
			case(5):y = 1; z = 3; break;
			case(6):y = 1; z = 2; break;
			}
			mat.cijkl_glb[i - 1][j - 1] = Cwxyz[w - 1][x - 1][y - 1][z - 1];
		}
	}
}

void global_parameters::copy_to_device() {
#pragma acc enter data copyin(this)

#pragma acc enter data copyin(this->Hext[0:3])
#pragma acc enter data copyin(this->Hext_stat[0:3])
#pragma acc enter data copyin(this->Hext_altr[0:3])
#pragma acc enter data copyin(this->H_altr_freq[0:3])

#pragma acc enter data copyin(this->Eext[0:3])
#pragma acc enter data copyin(this->Eext_stat[0:3])
#pragma acc enter data copyin(this->Eext_altr[0:3])
#pragma acc enter data copyin(this->E_altr_freq[0:3])

#pragma acc enter data copyin(this->planeEM_E[0:3])
#pragma acc enter data copyin(this->planeEM_E_freq[0:3])

#pragma acc enter data copyin(this->Tg2c[0:9])
#pragma acc enter data copyin(this->Tc2g[0:9])

#pragma acc enter data copyin(this->material_parameters[0:this->num_materials])

#pragma acc enter data copyin(this->material_cell)
#pragma acc enter data copyin(this->material_cell.matrix[0:this->n])

#pragma acc enter data copyin(this->free_charge_surfaces[0:nz+1])
#pragma acc enter data copyin(this->free_charge_cells[0:nz])

	if (if_output_only_magcell == true) {
#pragma acc enter data copyin(this->magcell_index[0:NFM])
	}

	if (if_spin_pumping == true) {
#pragma acc enter data copyin(this->if_J_ISHE[0:nz])
#pragma acc enter data copyin(this->source_z[0:nz])
#pragma acc enter data copyin(this->thickness_pumping_layer[0:nz])
	}

	if (if_prescribe_Eext == true) {
#pragma acc enter data copyin(this->prescribe_Eext_time[0:num_prescribe_Eext])
#pragma acc enter data copyin(this->prescribe_Ex[0:num_prescribe_Eext])
#pragma acc enter data copyin(this->prescribe_Ey[0:num_prescribe_Eext])
#pragma acc enter data copyin(this->prescribe_Ez[0:num_prescribe_Eext])
	}

	if (if_prescribe_Hext == true) {
#pragma acc enter data copyin(this->prescribe_Hext_time[0:num_prescribe_Hext])
#pragma acc enter data copyin(this->prescribe_Hx[0:num_prescribe_Hext])
#pragma acc enter data copyin(this->prescribe_Hy[0:num_prescribe_Hext])
#pragma acc enter data copyin(this->prescribe_Hz[0:num_prescribe_Hext])
	}

	{ //Non-uniform external magnetic field
#pragma acc enter data copyin(this->Hextx_nonunif)
#pragma acc enter data copyin(this->Hexty_nonunif)
#pragma acc enter data copyin(this->Hextz_nonunif)
#pragma acc enter data copyin(this->Hextx_nonunif.matrix[0:this->n])
#pragma acc enter data copyin(this->Hexty_nonunif.matrix[0:this->n])
#pragma acc enter data copyin(this->Hextz_nonunif.matrix[0:this->n])	
	}

#pragma acc enter data copyin(this->kx[0:this->nx], this->ky[0:this->ny], this->kz[0:(this->nz)/2+1])
#pragma acc enter data copyin(this->k_norm)
#pragma acc enter data copyin(this->k_norm.matrix[0:(this->nx)*(this->ny)*((this->nz)/2+1)])
}

void global_parameters::update_to_device() {
#pragma acc update device(this->Hext[0:3])
}

#pragma acc routine seq nohost
void global_parameters::time_marching_device_serial() {
	time_device = time_device + dt;

	if (if_prescribe_Eext == true) {
		if ((time_device - dt) > prescribe_Eext_time[prescribe_index_E]) {
			prescribe_index_E = prescribe_index_E + 1;
		}

		if (prescribe_index_E < num_prescribe_Eext) {
			Eext[0] = prescribe_Ex[prescribe_index_E];
			Eext[1] = prescribe_Ey[prescribe_index_E];
			Eext[2] = prescribe_Ez[prescribe_index_E];
		}
		else {
			Eext[0] = 0.;
			Eext[1] = 0.;
			Eext[2] = 0.;
		}
	}

	if (if_prescribe_Hext == true) {
		if ((time_device - dt) > prescribe_Hext_time[prescribe_index_H]) {
			prescribe_index_H = prescribe_index_H + 1;
		}

		if (prescribe_index_H < num_prescribe_Hext) {
			Hext[0] = prescribe_Hx[prescribe_index_H];
			Hext[1] = prescribe_Hy[prescribe_index_H];
			Hext[2] = prescribe_Hz[prescribe_index_H];
		}
		else {
			Hext[0] = 0.;
			Hext[1] = 0.;
			Hext[2] = 0.;
		}
	}
}

global_parameters global_parameters::glb;
