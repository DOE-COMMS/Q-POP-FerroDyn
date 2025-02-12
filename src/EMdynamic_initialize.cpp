#include "EMdynamic_system.h"

void EMdynamic_system::initialize_host(){

//	pt_mag=pt_mag_in;
//	pt_ele=pt_ele_in;

	nx = pt_geo->nx_system;
	ny = pt_geo->ny_system;
	nz = pt_geo->nz_system;

	n = nx * ny * nz;

	std::cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;

	dx = pt_geo->dx;
	dy = pt_geo->dy;
	dz = pt_geo->dz;
	
	nx_phy = pt_geo->nx_phy;
	ny_phy = pt_geo->ny_phy;
	nz_phy = pt_geo->nz_phy;

	xS = pt_geo->xS; yS = pt_geo->yS; zS = pt_geo->zS;
	xE = pt_geo->xE; yE = pt_geo->yE; zE = pt_geo->zE;

	if_PML = pt_geo->if_PML;
	if_PML_Xs = pt_geo->if_PML_Xs; if_PML_Xe = pt_geo->if_PML_Xe;
	if_PML_Ys = pt_geo->if_PML_Ys; if_PML_Ye = pt_geo->if_PML_Ye;
	if_PML_Zs = pt_geo->if_PML_Zs; if_PML_Ze = pt_geo->if_PML_Ze;

	periodicX_EM = pt_geo->periodicX_EM;
	periodicY_EM = pt_geo->periodicY_EM;
	periodicZ_EM = pt_geo->periodicZ_EM;

	PML_size = pt_geo->PML_size;

	planeEM_source_z = pt_geo->pt_nz_layer[0] / 5;

	int i, j;
	unsigned int mat_type;
	material* mat;

	sqrt_er_bot = 1.;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (pt_glb->material_cell(i, j, pt_geo->zS) != 0) {
				mat_type = pt_glb->material_cell(i, j, pt_geo->zS);
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				sqrt_er_bot = sqrt(mat->r_permittivity11);

				j = ny; i = nx;
			}
		}
	}

	sqrt_er_top = 1.;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (pt_glb->material_cell(i, j, pt_geo->zE - 1) != 0) {
				mat_type = pt_glb->material_cell(i, j, pt_geo->zE - 1);
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				sqrt_er_top = sqrt(mat->r_permittivity11);

				j = ny; i = nx;
			}
		}
	}

	if (pt_glb->if_Jf_input == true) {		
		Jf_nx = (pt_glb->Jfin_xf - pt_glb->Jfin_xi + 1);
		Jf_ny = (pt_glb->Jfin_yf - pt_glb->Jfin_yi + 1);
		Jf_nz = (pt_glb->Jfin_zf - pt_glb->Jfin_zi + 1);
		Jf_n = Jf_nx * Jf_ny * Jf_nz;
	}

	//Initialization
	DHx_em.initialize(nx + 1, ny, nz); DHx_em_store.initialize(nx + 1, ny, nz);
	DHy_em.initialize(nx, ny + 1, nz); DHy_em_store.initialize(nx, ny + 1, nz);
	DHz_em.initialize(nx, ny, nz + 1); DHz_em_store.initialize(nx, ny, nz + 1);

	dDHx_em_rk1.initialize(nx + 1, ny, nz);
	dDHy_em_rk1.initialize(nx, ny + 1, nz);
	dDHz_em_rk1.initialize(nx, ny, nz + 1);

	dDHx_em_rk2.initialize(nx + 1, ny, nz);
	dDHy_em_rk2.initialize(nx, ny + 1, nz);
	dDHz_em_rk2.initialize(nx, ny, nz + 1);

	dDHx_em_rk3.initialize(nx + 1, ny, nz);
	dDHy_em_rk3.initialize(nx, ny + 1, nz);
	dDHz_em_rk3.initialize(nx, ny, nz + 1);

	dDHx_em_rk4.initialize(nx + 1, ny, nz);
	dDHy_em_rk4.initialize(nx, ny + 1, nz);
	dDHz_em_rk4.initialize(nx, ny, nz + 1);

	DEx_em.initialize(nx, ny + 1, nz + 1); DEx_em_store.initialize(nx, ny + 1, nz + 1);
	DEy_em.initialize(nx + 1, ny, nz + 1); DEy_em_store.initialize(nx + 1, ny, nz + 1);
	DEz_em.initialize(nx + 1, ny + 1, nz); DEz_em_store.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk1.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk1.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk1.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk2.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk2.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk2.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk3.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk3.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk3.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk4.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk4.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk4.initialize(nx + 1, ny + 1, nz);

	DHx_em_cell.initialize(nx, ny, nz); DHy_em_cell.initialize(nx, ny, nz); DHz_em_cell.initialize(nx, ny, nz);
	DEx_em_cell.initialize(nx, ny, nz); DEy_em_cell.initialize(nx, ny, nz); DEz_em_cell.initialize(nx, ny, nz);

	DJfx.initialize(nx, ny, nz); DJfy.initialize(nx, ny, nz); DJfz.initialize(nx, ny, nz);

	Jpx_n1.initialize(nx, ny, nz); Jpy_n1.initialize(nx, ny, nz); Jpz_n1.initialize(nx, ny, nz);
	Jpx_n2.initialize(nx, ny, nz); Jpy_n2.initialize(nx, ny, nz); Jpz_n2.initialize(nx, ny, nz);
	Jpx_n3.initialize(nx, ny, nz); Jpy_n3.initialize(nx, ny, nz); Jpz_n3.initialize(nx, ny, nz);

	Jpx_n1_store.initialize(nx, ny, nz); Jpy_n1_store.initialize(nx, ny, nz); Jpz_n1_store.initialize(nx, ny, nz);
	Jpx_n2_store.initialize(nx, ny, nz); Jpy_n2_store.initialize(nx, ny, nz); Jpz_n2_store.initialize(nx, ny, nz);
	Jpx_n3_store.initialize(nx, ny, nz); Jpy_n3_store.initialize(nx, ny, nz); Jpz_n3_store.initialize(nx, ny, nz);

	dJpx_n1_rk1.initialize(nx, ny, nz); dJpy_n1_rk1.initialize(nx, ny, nz); dJpz_n1_rk1.initialize(nx, ny, nz);
	dJpx_n2_rk1.initialize(nx, ny, nz); dJpy_n2_rk1.initialize(nx, ny, nz); dJpz_n2_rk1.initialize(nx, ny, nz);
	dJpx_n3_rk1.initialize(nx, ny, nz); dJpy_n3_rk1.initialize(nx, ny, nz); dJpz_n3_rk1.initialize(nx, ny, nz);

	dJpx_n1_rk2.initialize(nx, ny, nz); dJpy_n1_rk2.initialize(nx, ny, nz); dJpz_n1_rk2.initialize(nx, ny, nz);
	dJpx_n2_rk2.initialize(nx, ny, nz); dJpy_n2_rk2.initialize(nx, ny, nz); dJpz_n2_rk2.initialize(nx, ny, nz);
	dJpx_n3_rk2.initialize(nx, ny, nz); dJpy_n3_rk2.initialize(nx, ny, nz); dJpz_n3_rk2.initialize(nx, ny, nz);

	dJpx_n1_rk3.initialize(nx, ny, nz); dJpy_n1_rk3.initialize(nx, ny, nz); dJpz_n1_rk3.initialize(nx, ny, nz);
	dJpx_n2_rk3.initialize(nx, ny, nz); dJpy_n2_rk3.initialize(nx, ny, nz); dJpz_n2_rk3.initialize(nx, ny, nz);
	dJpx_n3_rk3.initialize(nx, ny, nz); dJpy_n3_rk3.initialize(nx, ny, nz); dJpz_n3_rk3.initialize(nx, ny, nz);

	dJpx_n1_rk4.initialize(nx, ny, nz); dJpy_n1_rk4.initialize(nx, ny, nz); dJpz_n1_rk4.initialize(nx, ny, nz);
	dJpx_n2_rk4.initialize(nx, ny, nz); dJpy_n2_rk4.initialize(nx, ny, nz); dJpz_n2_rk4.initialize(nx, ny, nz);
	dJpx_n3_rk4.initialize(nx, ny, nz); dJpy_n3_rk4.initialize(nx, ny, nz); dJpz_n3_rk4.initialize(nx, ny, nz);

	DEx_em_t1.initialize(nx, ny + 1, nz + 1); DEy_em_t1.initialize(nx + 1, ny, nz + 1); DEz_em_t1.initialize(nx + 1, ny + 1, nz);
	DEx_em_t2.initialize(nx, ny + 1, nz + 1); DEy_em_t2.initialize(nx + 1, ny, nz + 1); DEz_em_t2.initialize(nx + 1, ny + 1, nz);
	DEx_em_t3.initialize(nx, ny + 1, nz + 1); DEy_em_t3.initialize(nx + 1, ny, nz + 1); DEz_em_t3.initialize(nx + 1, ny + 1, nz);
	DEx_em_t4.initialize(nx, ny + 1, nz + 1); DEy_em_t4.initialize(nx + 1, ny, nz + 1); DEz_em_t4.initialize(nx + 1, ny + 1, nz);
	
	if (pt_geo->if_PML == true) {
		Dx_PML.initialize(nx, ny + 1, nz + 1); Dy_PML.initialize(nx + 1, ny, nz + 1); Dz_PML.initialize(nx + 1, ny + 1, nz);
		Dx_PML_store.initialize(nx, ny + 1, nz + 1); Dy_PML_store.initialize(nx + 1, ny, nz + 1); Dz_PML_store.initialize(nx + 1, ny + 1, nz);

		Bx_PML.initialize(nx + 1, ny, nz); By_PML.initialize(nx, ny + 1, nz); Bz_PML.initialize(nx, ny, nz + 1);
		Bx_PML_store.initialize(nx + 1, ny, nz); By_PML_store.initialize(nx, ny + 1, nz); Bz_PML_store.initialize(nx, ny, nz + 1);

		BHx_em_t1.initialize(nx + 1, ny, nz); BHy_em_t1.initialize(nx, ny + 1, nz); BHz_em_t1.initialize(nx, ny, nz + 1);
		BHx_em_t2.initialize(nx + 1, ny, nz); BHy_em_t2.initialize(nx, ny + 1, nz); BHz_em_t2.initialize(nx, ny, nz + 1);
		BHx_em_t3.initialize(nx + 1, ny, nz); BHy_em_t3.initialize(nx, ny + 1, nz); BHz_em_t3.initialize(nx, ny, nz + 1);
		BHx_em_t4.initialize(nx + 1, ny, nz); BHy_em_t4.initialize(nx, ny + 1, nz); BHz_em_t4.initialize(nx, ny, nz + 1);
	}
	
	sigma_x_n.initialize(nx, 1, 1); sigma_y_n.initialize(1, ny, 1); sigma_z_n.initialize(1, 1, nz);
	kappa_x_n.initialize(nx, 1, 1); kappa_y_n.initialize(1, ny, 1); kappa_z_n.initialize(1, 1, nz);

	sigma_x_np1.initialize(nx + 1, 1, 1); sigma_y_np1.initialize(1, ny + 1, 1); sigma_z_np1.initialize(1, 1, nz + 1);
	kappa_x_np1.initialize(nx+1, 1, 1); kappa_y_np1.initialize(1, ny+1, 1); kappa_z_np1.initialize(1, 1, nz+1);

	initialize_PML();
}

void EMdynamic_system::initialize_device() {
	if (pt_glb->if_input_em == true) {

#pragma acc parallel default(present) async(3)
		{
			DHx_em_store = DHx_em;
			DHy_em_store = DHy_em;
			DHz_em_store = DHz_em;
			DEx_em_store = DEx_em;
			DEy_em_store = DEy_em;
			DEz_em_store = DEz_em;

			Dx_PML = DEx_em;
			Dy_PML = DEy_em;
			Dz_PML = DEz_em;
			Bx_PML = DHx_em;
			By_PML = DHy_em;
			Bz_PML = DHz_em;	

			Dx_PML_store = DEx_em;
			Dy_PML_store = DEy_em;
			Dz_PML_store = DEz_em;
			Bx_PML_store = DHx_em;
			By_PML_store = DHy_em;
			Bz_PML_store = DHz_em;
		}

#pragma acc parallel default(present) async(3)
		{
			update_DH_cell();
			update_DE_cell();
		}
	}

	if (pt_glb->if_input_Jp == true) {
#pragma acc parallel default(present) async(4)
		{
			Jpx_n1_store = Jpx_n1;
			Jpy_n1_store = Jpy_n1;
			Jpz_n1_store = Jpz_n1;
			Jpx_n2_store = Jpx_n2;
			Jpy_n2_store = Jpy_n2;
			Jpz_n2_store = Jpz_n2;
			Jpx_n3_store = Jpx_n3;
			Jpy_n3_store = Jpy_n3;
			Jpz_n3_store = Jpz_n3;
		}
	}
}

void EMdynamic_system::copy_to_device() {
#pragma acc enter data copyin(this)

#pragma acc serial default(present)
	{
		this->pt_geo = &(geometry_parameters::geo);
		this->pt_glb = &(global_parameters::glb);
		this->pt_math = &(mathlib::mlb);
	}

#pragma acc enter data copyin(this->DHx_em, this->DHx_em_store)
#pragma acc enter data copyin(this->DHy_em, this->DHy_em_store)
#pragma acc enter data copyin(this->DHz_em, this->DHz_em_store)

#pragma acc enter data copyin(this->DHx_em.matrix[0:(nx+1)*ny*nz], this->DHx_em_store.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->DHy_em.matrix[0:nx*(ny+1)*nz], this->DHy_em_store.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->DHz_em.matrix[0:nx*ny*(nz+1)], this->DHz_em_store.matrix[0:nx*ny*(nz+1)])

#pragma acc enter data copyin(this->DEx_em, this->DEx_em_store)
#pragma acc enter data copyin(this->DEy_em, this->DEy_em_store)
#pragma acc enter data copyin(this->DEz_em, this->DEz_em_store)

#pragma acc enter data copyin(this->DEx_em.matrix[0:nx*(ny+1)*(nz+1)], this->DEx_em_store.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em.matrix[0:(nx+1)*ny*(nz+1)], this->DEy_em_store.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em.matrix[0:(nx+1)*(ny+1)*nz], this->DEz_em_store.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->dDHx_em_rk1)
#pragma acc enter data copyin(this->dDHy_em_rk1)
#pragma acc enter data copyin(this->dDHz_em_rk1)
#pragma acc enter data copyin(this->dDHx_em_rk1.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk1.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk1.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk1)
#pragma acc enter data copyin(this->dDEy_em_rk1)
#pragma acc enter data copyin(this->dDEz_em_rk1)
#pragma acc enter data copyin(this->dDEx_em_rk1.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk1.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk1.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->dDHx_em_rk2)
#pragma acc enter data copyin(this->dDHy_em_rk2)
#pragma acc enter data copyin(this->dDHz_em_rk2)
#pragma acc enter data copyin(this->dDHx_em_rk2.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk2.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk2.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk2)
#pragma acc enter data copyin(this->dDEy_em_rk2)
#pragma acc enter data copyin(this->dDEz_em_rk2)
#pragma acc enter data copyin(this->dDEx_em_rk2.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk2.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk2.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->dDHx_em_rk3)
#pragma acc enter data copyin(this->dDHy_em_rk3)
#pragma acc enter data copyin(this->dDHz_em_rk3)
#pragma acc enter data copyin(this->dDHx_em_rk3.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk3.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk3.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk3)
#pragma acc enter data copyin(this->dDEy_em_rk3)
#pragma acc enter data copyin(this->dDEz_em_rk3)
#pragma acc enter data copyin(this->dDEx_em_rk3.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk3.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk3.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->dDHx_em_rk4)
#pragma acc enter data copyin(this->dDHy_em_rk4)
#pragma acc enter data copyin(this->dDHz_em_rk4)
#pragma acc enter data copyin(this->dDHx_em_rk4.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk4.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk4.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk4)
#pragma acc enter data copyin(this->dDEy_em_rk4)
#pragma acc enter data copyin(this->dDEz_em_rk4)
#pragma acc enter data copyin(this->dDEx_em_rk4.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk4.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk4.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->DHx_em_cell,this->DHy_em_cell,this->DHz_em_cell)
#pragma acc enter data copyin(this->DEx_em_cell,this->DEy_em_cell,this->DEz_em_cell)

#pragma acc enter data copyin(this->DHx_em_cell.matrix[0:n],this->DHy_em_cell.matrix[0:n],this->DHz_em_cell.matrix[0:n])
#pragma acc enter data copyin(this->DEx_em_cell.matrix[0:n],this->DEy_em_cell.matrix[0:n],this->DEz_em_cell.matrix[0:n])

#pragma acc enter data copyin(this->DJfx,this->DJfy,this->DJfz)
#pragma acc enter data copyin(this->DJfx.matrix[0:n],this->DJfy.matrix[0:n],this->DJfz.matrix[0:n])

#pragma acc enter data copyin(this->Jpx_n1,this->Jpy_n1,this->Jpz_n1)
#pragma acc enter data copyin(this->Jpx_n2,this->Jpy_n2,this->Jpz_n2)
#pragma acc enter data copyin(this->Jpx_n3,this->Jpy_n3,this->Jpz_n3)
#pragma acc enter data copyin(this->Jpx_n1_store,this->Jpy_n1_store,this->Jpz_n1_store)
#pragma acc enter data copyin(this->Jpx_n2_store,this->Jpy_n2_store,this->Jpz_n2_store)
#pragma acc enter data copyin(this->Jpx_n3_store,this->Jpy_n3_store,this->Jpz_n3_store)

#pragma acc enter data copyin(this->Jpx_n1.matrix[0:n],this->Jpy_n1.matrix[0:n],this->Jpz_n1.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n2.matrix[0:n],this->Jpy_n2.matrix[0:n],this->Jpz_n2.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n3.matrix[0:n],this->Jpy_n3.matrix[0:n],this->Jpz_n3.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n1_store.matrix[0:n],this->Jpy_n1_store.matrix[0:n],this->Jpz_n1_store.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n2_store.matrix[0:n],this->Jpy_n2_store.matrix[0:n],this->Jpz_n2_store.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n3_store.matrix[0:n],this->Jpy_n3_store.matrix[0:n],this->Jpz_n3_store.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk1,			 this->dJpy_n1_rk1,				this->dJpz_n1_rk1)
#pragma acc enter data copyin(this->dJpx_n2_rk1,			 this->dJpy_n2_rk1,				this->dJpz_n2_rk1)
#pragma acc enter data copyin(this->dJpx_n3_rk1,			 this->dJpy_n3_rk1,				this->dJpz_n3_rk1)
#pragma acc enter data copyin(this->dJpx_n1_rk1.matrix[0:n], this->dJpy_n1_rk1.matrix[0:n], this->dJpz_n1_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk1.matrix[0:n], this->dJpy_n2_rk1.matrix[0:n], this->dJpz_n2_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk1.matrix[0:n], this->dJpy_n3_rk1.matrix[0:n], this->dJpz_n3_rk1.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk2,			 this->dJpy_n1_rk2,				this->dJpz_n1_rk2)
#pragma acc enter data copyin(this->dJpx_n2_rk2,			 this->dJpy_n2_rk2,				this->dJpz_n2_rk2)
#pragma acc enter data copyin(this->dJpx_n3_rk2,			 this->dJpy_n3_rk2,				this->dJpz_n3_rk2)
#pragma acc enter data copyin(this->dJpx_n1_rk2.matrix[0:n], this->dJpy_n1_rk2.matrix[0:n], this->dJpz_n1_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk2.matrix[0:n], this->dJpy_n2_rk2.matrix[0:n], this->dJpz_n2_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk2.matrix[0:n], this->dJpy_n3_rk2.matrix[0:n], this->dJpz_n3_rk2.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk3,			 this->dJpy_n1_rk3,				this->dJpz_n1_rk3)
#pragma acc enter data copyin(this->dJpx_n2_rk3,			 this->dJpy_n2_rk3,				this->dJpz_n2_rk3)
#pragma acc enter data copyin(this->dJpx_n3_rk3,			 this->dJpy_n3_rk3,				this->dJpz_n3_rk3)
#pragma acc enter data copyin(this->dJpx_n1_rk3.matrix[0:n], this->dJpy_n1_rk3.matrix[0:n], this->dJpz_n1_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk3.matrix[0:n], this->dJpy_n2_rk3.matrix[0:n], this->dJpz_n2_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk3.matrix[0:n], this->dJpy_n3_rk3.matrix[0:n], this->dJpz_n3_rk3.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk4,			 this->dJpy_n1_rk4,				this->dJpz_n1_rk4)
#pragma acc enter data copyin(this->dJpx_n2_rk4,			 this->dJpy_n2_rk4,				this->dJpz_n2_rk4)
#pragma acc enter data copyin(this->dJpx_n3_rk4,			 this->dJpy_n3_rk4,				this->dJpz_n3_rk4)
#pragma acc enter data copyin(this->dJpx_n1_rk4.matrix[0:n], this->dJpy_n1_rk4.matrix[0:n], this->dJpz_n1_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk4.matrix[0:n], this->dJpy_n2_rk4.matrix[0:n], this->dJpz_n2_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk4.matrix[0:n], this->dJpy_n3_rk4.matrix[0:n], this->dJpz_n3_rk4.matrix[0:n])

#pragma acc enter data copyin(this->DEx_em_t1)
#pragma acc enter data copyin(this->DEy_em_t1)
#pragma acc enter data copyin(this->DEz_em_t1)
			
#pragma acc enter data copyin(this->DEx_em_t2)
#pragma acc enter data copyin(this->DEy_em_t2)
#pragma acc enter data copyin(this->DEz_em_t2)

#pragma acc enter data copyin(this->DEx_em_t3)
#pragma acc enter data copyin(this->DEy_em_t3)
#pragma acc enter data copyin(this->DEz_em_t3)

#pragma acc enter data copyin(this->DEx_em_t4)
#pragma acc enter data copyin(this->DEy_em_t4)
#pragma acc enter data copyin(this->DEz_em_t4)


#pragma acc enter data copyin(this->DEx_em_t1.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t1.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t1.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->DEx_em_t2.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t2.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t2.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->DEx_em_t3.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t3.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t3.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->DEx_em_t4.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t4.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t4.matrix[0:(nx+1)*(ny+1)*nz])

if (pt_geo->if_PML == true) {
		
#pragma acc enter data copyin(this->Dx_PML, this->Dy_PML, this->Dz_PML)
#pragma acc enter data copyin(this->Dx_PML_store, this->Dy_PML_store, this->Dz_PML_store)

#pragma acc enter data copyin(this->Bx_PML, this->By_PML, this->Bz_PML)
#pragma acc enter data copyin(this->Bx_PML_store, this->By_PML_store, this->Bz_PML_store)
			
#pragma acc enter data copyin(this->BHx_em_t1)
#pragma acc enter data copyin(this->BHy_em_t1)
#pragma acc enter data copyin(this->BHz_em_t1)

#pragma acc enter data copyin(this->BHx_em_t2)
#pragma acc enter data copyin(this->BHy_em_t2)
#pragma acc enter data copyin(this->BHz_em_t2)

#pragma acc enter data copyin(this->BHx_em_t3)
#pragma acc enter data copyin(this->BHy_em_t3)
#pragma acc enter data copyin(this->BHz_em_t3)

#pragma acc enter data copyin(this->BHx_em_t4)
#pragma acc enter data copyin(this->BHy_em_t4)
#pragma acc enter data copyin(this->BHz_em_t4)

#pragma acc enter data copyin(this->Dx_PML.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->Dy_PML.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->Dz_PML.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->Dx_PML_store.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->Dy_PML_store.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->Dz_PML_store.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->Bx_PML.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->By_PML.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->Bz_PML.matrix[0:nx*ny*(nz+1)])

#pragma acc enter data copyin(this->Bx_PML_store.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->By_PML_store.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->Bz_PML_store.matrix[0:nx*ny*(nz+1)])

#pragma acc enter data copyin(this->BHx_em_t1.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->BHy_em_t1.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->BHz_em_t1.matrix[0:nx*ny*(nz+1)])

#pragma acc enter data copyin(this->BHx_em_t2.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->BHy_em_t2.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->BHz_em_t2.matrix[0:nx*ny*(nz+1)])

#pragma acc enter data copyin(this->BHx_em_t3.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->BHy_em_t3.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->BHz_em_t3.matrix[0:nx*ny*(nz+1)])

#pragma acc enter data copyin(this->BHx_em_t4.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->BHy_em_t4.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->BHz_em_t4.matrix[0:nx*ny*(nz+1)])	

#pragma acc enter data copyin(this->sigma_x_n)
#pragma acc enter data copyin(this->sigma_y_n)
#pragma acc enter data copyin(this->sigma_z_n)
#pragma acc enter data copyin(this->sigma_x_n.matrix[0:nx])
#pragma acc enter data copyin(this->sigma_y_n.matrix[0:ny])
#pragma acc enter data copyin(this->sigma_z_n.matrix[0:nz])

#pragma acc enter data copyin(this->kappa_x_n)
#pragma acc enter data copyin(this->kappa_y_n)
#pragma acc enter data copyin(this->kappa_z_n)
#pragma acc enter data copyin(this->kappa_x_n.matrix[0:nx])
#pragma acc enter data copyin(this->kappa_y_n.matrix[0:ny])
#pragma acc enter data copyin(this->kappa_z_n.matrix[0:nz])

#pragma acc enter data copyin(this->sigma_x_np1)
#pragma acc enter data copyin(this->sigma_y_np1)
#pragma acc enter data copyin(this->sigma_z_np1)
#pragma acc enter data copyin(this->sigma_x_np1.matrix[0:nx+1])
#pragma acc enter data copyin(this->sigma_y_np1.matrix[0:ny+1])
#pragma acc enter data copyin(this->sigma_z_np1.matrix[0:nz+1])


#pragma acc enter data copyin(this->kappa_x_np1)
#pragma acc enter data copyin(this->kappa_y_np1)
#pragma acc enter data copyin(this->kappa_z_np1)
#pragma acc enter data copyin(this->kappa_x_np1.matrix[0:nx+1])
#pragma acc enter data copyin(this->kappa_y_np1.matrix[0:ny+1])
#pragma acc enter data copyin(this->kappa_z_np1.matrix[0:nz+1])
		}
}

void EMdynamic_system::copy_from_device() {
#pragma acc update host(DHx_em_cell.matrix[0:n],DHy_em_cell.matrix[0:n],DHz_em_cell.matrix[0:n])
#pragma acc update host(DEx_em_cell.matrix[0:n],DEy_em_cell.matrix[0:n],DEz_em_cell.matrix[0:n])
}

void EMdynamic_system::copyYee_from_device() {
#pragma acc update host(DHx_em.matrix[0:(nx+1)*ny*nz])
#pragma acc update host(DHy_em.matrix[0:nx*(ny+1)*nz])
#pragma acc update host(DHz_em.matrix[0:nx*ny*(nz+1)])

#pragma acc update host(DEx_em.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc update host(DEy_em.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc update host(DEz_em.matrix[0:(nx+1)*(ny+1)*nz])
}

void EMdynamic_system::copyJp_from_device() {
#pragma acc update host(Jpx_n1.matrix[0:n],Jpy_n1.matrix[0:n],Jpz_n1.matrix[0:n])
#pragma acc update host(Jpx_n2.matrix[0:n],Jpy_n2.matrix[0:n],Jpz_n2.matrix[0:n])
#pragma acc update host(Jpx_n3.matrix[0:n],Jpy_n3.matrix[0:n],Jpz_n3.matrix[0:n])
}

// Calculate PML parameters (sigma and kappa) at each point in the mesh
void EMdynamic_system::initialize_PML() {
	for (long int i = 0; i < nx; i++)
	{
		sigma_x_n(i) = 0.0;
		kappa_x_n(i) = 1.0;

		if (i < pt_geo->PML_size && pt_geo->if_PML_Xs)
		{
			sigma_x_n(i) = pow((double)(pt_geo->PML_size - i) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_x_n(i) = 1.0 + pow((double)(pt_geo->PML_size - i) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}
		if (i >= nx - pt_geo->PML_size && pt_geo->if_PML_Xe)
		{
			sigma_x_n(i) = pow((double)(pt_geo->PML_size - (nx - i - 1)) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_x_n(i) = 1.0 + pow((pt_geo->PML_size - (nx - i - 1)) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}

		//std::cout << "sigma_x_n(" << i << "): " << sigma_x_n(i) << std::endl;
		//std::cout << "kappa_x_n(" << i << "): " << kappa_x_n(i) << std::endl;
	}

	for (long int j = 0; j < ny; j++)
	{
		sigma_y_n(j) = 0.0;
		kappa_y_n(j) = 1.0;

		if (j < pt_geo->PML_size && pt_geo->if_PML_Ys)
		{
			sigma_y_n(j) = pow((double)(pt_geo->PML_size - j) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_y_n(j) = 1.0 + pow((double)(pt_geo->PML_size - j) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}
		if (j >= ny - pt_geo->PML_size && pt_geo->if_PML_Ye)
		{
			sigma_y_n(j) = pow((double)(pt_geo->PML_size - (ny - j - 1)) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_y_n(j) = 1.0 + pow((double)(pt_geo->PML_size - (ny - j - 1)) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}

		//std::cout << "sigma_y_n(" << j << "): " << sigma_y_n(j) << std::endl;
		//std::cout << "kappa_y_n(" << j << "): " << kappa_y_n(j) << std::endl;
	}

	for (long int k = 0; k < nz; k++)
	{
		sigma_z_n(k) = 0.0;
		kappa_z_n(k) = 1.0;

		if (k < pt_geo->PML_size && pt_geo->if_PML_Zs)
		{
			sigma_z_n(k) = pow((double)(pt_geo->PML_size - k) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_z_n(k) = 1.0 + pow((double)(pt_geo->PML_size - k) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}
		if (k >= nz - pt_geo->PML_size && pt_geo->if_PML_Ze)
		{
			sigma_z_n(k) = pow((double)(pt_geo->PML_size - (nz - k - 1)) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_z_n(k) = 1.0 + pow((double)(pt_geo->PML_size - (nz - k - 1)) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}

		//std::cout << "sigma_z_n(" << k << "): " << sigma_z_n(k) << std::endl;
		//std::cout << "kappa_z_n(" << k << "): " << kappa_z_n(k) << std::endl;
	}

	for (long int i = 0; i < nx+1; i++)
	{
		sigma_x_np1(i) = 0.0;
		kappa_x_np1(i) = 1.0;

		if (i < pt_geo->PML_size && pt_geo->if_PML_Xs)
		{
			sigma_x_np1(i) = pow(((double)(pt_geo->PML_size - i) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * \
				pt_glb->sigmaMax;
			kappa_x_np1(i) = 1.0 + pow(((double)(pt_geo->PML_size - i) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}
		if (i >= nx + 1 - pt_geo->PML_size && pt_geo->if_PML_Xe)
		{
			sigma_x_np1(i) = pow(((double)(pt_geo->PML_size - (nx - i)) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * \
				pt_glb->sigmaMax;
			kappa_x_np1(i) = 1.0 + pow(((double)(pt_geo->PML_size - (nx - i)) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}

		//std::cout << "sigma_x_np1(" << i << "): " << sigma_x_np1(i) << std::endl;
		//std::cout << "kappa_x_np1(" << i << "): " << kappa_x_np1(i) << std::endl;
	}

	for (long int j = 0; j < ny + 1; j++)
	{
		sigma_y_np1(j) = 0.0;
		kappa_y_np1(j) = 1.0;

		if (j < pt_geo->PML_size && pt_geo->if_PML_Ys)
		{
			sigma_y_np1(j) = pow(((double)(pt_geo->PML_size - j) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_y_np1(j) = 1.0 + pow(((double)(pt_geo->PML_size - j) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}
		if (j >= ny + 1 - pt_geo->PML_size && pt_geo->if_PML_Ye)
		{
			sigma_y_np1(j) = pow(((double)(pt_geo->PML_size - (ny - j)) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_y_np1(j) = 1.0 + pow(((double)(pt_geo->PML_size - (ny - j)) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}

		//std::cout << "sigma_y_np1(" << j << "): " << sigma_y_np1(j) << std::endl;
		//std::cout << "kappa_y_np1(" << j << "): " << kappa_y_np1(j) << std::endl;
	}

	for (long int k = 0; k < nz + 1; k++)
	{
		sigma_z_np1(k) = 0.0;
		kappa_z_np1(k) = 1.0;

		if (k < pt_geo->PML_size && pt_geo->if_PML_Zs)
		{
			sigma_z_np1(k) = pow(((double)(pt_geo->PML_size - k) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_z_np1(k) = 1.0 + pow(((double)(pt_geo->PML_size - k) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}
		if (k >= nz + 1 - pt_geo->PML_size && pt_geo->if_PML_Ze)
		{
			sigma_z_np1(k) = pow(((double)(pt_geo->PML_size - (nz - k)) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) \
				* pt_glb->sigmaMax;
			kappa_z_np1(k) = 1.0 + pow(((double)(pt_geo->PML_size - (nz - k)) + 0.5) / (double)pt_geo->PML_size, pt_glb->PML_m) * (pt_glb->kappaMax - 1.0);
		}

		//std::cout << "sigma_z_np1(" << k << "): " << sigma_z_np1(k) << std::endl;
		//std::cout << "kappa_z_np1(" << k << "): " << kappa_z_np1(k) << std::endl;
	}
}
