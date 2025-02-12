#include "EMdynamic_system.h"
#define CHECK_AND_PRINT_NAN(var, i, j, k) if (isnan(var)) { \
    printf("NaN detected in %s at (%d, %d, %d)\n", #var, i, j, k); \
    }

void EMdynamic_system::get_dE_RK1() {

	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;

	/*
	* 1: (i-1, j-1, k-1)
	* 2: (i-1, j-1, k)
	* 3: (i-1, j, k-1)
	* 4: (i-1, j, k)
	* 5: (i, j-1, k-1)
	* 6: (i, j-1, k)
	* 7: (i, j, k-1)
	* 8: (i, j, k)
	*/
	long int idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8;
	long int idy1, idy2, idy3, idy4, idy5, idy6, idy7, idy8;
	long int idz1, idz2, idz3, idz4, idz5, idz6, idz7, idz8;

	double er11, er12, er13, er21, er22, er23, er31, er32, er33;

	double Denominator;
	double inverse_er11, inverse_er12, inverse_er13;
	double inverse_er21, inverse_er22, inverse_er23;
	double inverse_er31, inverse_er32, inverse_er33;

	double cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	bool isPML = false;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK1();
	}

#pragma acc wait


	//-----------X Y Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1, idx2, idy2, idz2, idx3, idy3, idz3, idx4, idy4, idz4,\
idx5, idy5, idz5, idx6, idy6, idz6, idx7, idy7, idz7, idx8, idy8, idz8,\
er11, er12, er13, er21, er22, er23, er31, er32, er33,\
cond11, con12, cond13, cond21, con22, cond23, cond31, con32, cond33,\
Denominator,\
inverse_er11, inverse_er12, inverse_er13,\
inverse_er21, inverse_er22, inverse_er23,\
inverse_er31, inverse_er32, inverse_er33,\
cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,isPML,  \
P_count) default(present)async(1)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

		if ((i == 0 || i == nx) && periodicX_EM == false) {
			if (if_PML_Xe == false && if_PML_Xs == false) {
				continue;
			}
			else {
				idx1 = i; idx2 = i; idx3 = i; idx4 = i; idx5 = i; idx6 = i; idx7 = i; idx8 = i;	
			}
		}
		else if ((i == 0 || i == nx) && periodicX_EM == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = nx - 1; idx4 = nx - 1; idx5 = 0; idx6 = 0; idx7 = 0; idx8 = 0;
		}
		else if (i < nx + 1) {
			idx1 = i - 1; idx2 = i - 1; idx3 = i - 1; idx4 = i - 1; idx5 = i; idx6 = i; idx7 = i; idx8 = i;
		}

		if ((j == 0 || j == ny) && periodicY_EM == false) {
			if (if_PML_Ye == false && if_PML_Ys == false) {
				continue;
			}
			else {
			idy1 = j; idy2 = j; idy3 = j; idy4 = j; idy5 = j; idy6 = j; idy7 = j; idy8 = j;
			}
		}
		else if ((j == 0 || j == ny) && periodicY_EM == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0; idy5 = ny - 1; idy6 = ny - 1; idy7 = 0; idy8 = 0;
		}
		else if (j < ny + 1) {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j; idy5 = j - 1; idy6 = j - 1; idy7 = j; idy8 = j;
		}
		
		if ((k == 0 || k == nz) && periodicZ_EM == false) {
			if (if_PML_Ze == false && if_PML_Zs == false) {
				continue;
			}
			else {
				idz1 = k; idz2 = k; idz3 = k; idz4 = k; idz5 = k; idz6 = k; idz7 = k; idz8 = k;
			}
		}
		else if ((k == 0 || k == nz) && periodicZ_EM == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0; idz5 = nz - 1; idz6 = 0; idz7 = nz - 1; idz8 = 0;
		}
		else if (k < nz + 1) {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k; idz5 = k - 1; idz6 = k; idz7 = k - 1; idz8 = k;
		}

		dHzdy = (DHz_em_store(idx5, idy3, idz4) - DHz_em_store(idx5, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx5, idy3, idz1)) / dz;
		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.;

		dHxdz = (DHx_em_store(idx5, idy7, idz2) - DHx_em_store(idx5, idy7, idz1)) / dz;
		dHzdx = (DHz_em_store(idx5, idy7, idz2) - DHz_em_store(idx1, idy7, idz2)) / dx;
		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.;

		dHydx = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx1, idy3, idz2)) / dx;
		dHxdy = (DHx_em_store(idx5, idy3, idz2) - DHx_em_store(idx5, idy1, idz2)) / dy;
		Jfz = 0.; Jpz = 0.;
		dPz = 0.;

		cond11 = 0.; cond12 = 0.; cond13 = 0.; cond21 = 0.; cond22 = 0.; cond23 = 0.; cond31 = 0.; cond32 = 0.; cond33 = 0.;
		er11 = 0.; er12 = 0.; er13 = 0.; er21 = 0.; er22 = 0.; er23 = 0.; er31 = 0.; er32 = 0.; er33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;

		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			if ((idx1 < xS && true == if_PML_Xs) \
			|| (idx1 >= xE && true == if_PML_Xe) \
			|| (idy1 < yS && true == if_PML_Ys) \
			|| (idy1 >= yE && true == if_PML_Ye) \
			|| (idz1 < zS && true == if_PML_Zs) \
			|| (idz1 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx1, idy1, idz1); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx1, idy1, idz1); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}

			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx1, idy1, idz1);
		Jfy = Jfy + DJfy(idx1, idy1, idz1);
		Jfz = Jfz + DJfz(idx1, idy1, idz1);
		Jpx = Jpx + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jpy = Jpy + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jpz = Jpz + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx1, idy1, idz1);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx1, idy1, idz1);


		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			if ((idx2 < xS && true == if_PML_Xs) \
				|| (idx2 >= xE && true == if_PML_Xe) \
				|| (idy2 < yS && true == if_PML_Ys) \
				|| (idy2 >= yE && true == if_PML_Ye) \
				|| (idz2 < zS && true == if_PML_Zs) \
				|| (idz2 >= zE && true == if_PML_Ze))
			{
				isPML = true;
			}
			else
				isPML = false;

			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx2, idy2, idz2); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx2, idy2, idz2); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx2, idy2, idz2);
		Jfy = Jfy + DJfy(idx2, idy2, idz2);
		Jfz = Jfz + DJfz(idx2, idy2, idz2);
		Jpx = Jpx + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jpy = Jpy + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jpz = Jpz + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx2, idy2, idz2);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx2, idy2, idz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		////printf("mat_type3 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx3, idy3, idz3);
		if (mat_type == 0) {
			if ((idx3 < xS && true == if_PML_Xs) \
			|| (idx3 >= xE && true == if_PML_Xe) \
			|| (idy3 < yS && true == if_PML_Ys) \
			|| (idy3 >= yE && true == if_PML_Ye) \
			|| (idz3 < zS && true == if_PML_Zs) \
			|| (idz3 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx3, idy3, idz3); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx3, idy3, idz3); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx3, idy3, idz3);
		Jfy = Jfy + DJfy(idx3, idy3, idz3);
		Jfz = Jfz + DJfz(idx3, idy3, idz3);
		Jpx = Jpx + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jpy = Jpy + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jpz = Jpz + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx3, idy3, idz3);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx3, idy3, idz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		//printf("mat_type4 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx4, idy4, idz4);
		if (mat_type == 0) {
			if ((idx4 < xS && true == if_PML_Xs) \
			|| (idx4 >= xE && true == if_PML_Xe) \
			|| (idy4 < yS && true == if_PML_Ys) \
			|| (idy4 >= yE && true == if_PML_Ye) \
			|| (idz4 < zS && true == if_PML_Zs) \
			|| (idz4 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx4, idy4, idz4); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx4, idy4, idz4); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx4, idy4, idz4);
		Jfy = Jfy + DJfy(idx4, idy4, idz4);
		Jfz = Jfz + DJfz(idx4, idy4, idz4);
		Jpx = Jpx + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jpy = Jpy + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jpz = Jpz + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx4, idy4, idz4);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx4, idy4, idz4);


		//----------Polarization 5------//
		mat_type = pt_glb->material_cell(idx5, idy5, idz5);
		if (mat_type == 0) {
			if ((idx5 < xS && true == if_PML_Xs) \
			|| (idx5 >= xE && true == if_PML_Xe) \
			|| (idy5 < yS && true == if_PML_Ys) \
			|| (idy5 >= yE && true == if_PML_Ye) \
			|| (idz5 < zS && true == if_PML_Zs) \
			|| (idz5 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx5, idy5, idz5); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx5, idy5, idz5); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx5, idy5, idz5); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx5, idy5, idz5);
		Jfy = Jfy + DJfy(idx5, idy5, idz5);
		Jfz = Jfz + DJfz(idx5, idy5, idz5);
		Jpx = Jpx + Jpx_n1_store(idx5, idy5, idz5) + Jpx_n2_store(idx5, idy5, idz5) + Jpx_n3_store(idx5, idy5, idz5);
		Jpy = Jpy + Jpy_n1_store(idx5, idy5, idz5) + Jpy_n2_store(idx5, idy5, idz5) + Jpy_n3_store(idx5, idy5, idz5);
		Jpz = Jpz + Jpz_n1_store(idx5, idy5, idz5) + Jpz_n2_store(idx5, idy5, idz5) + Jpz_n3_store(idx5, idy5, idz5);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx5, idy5, idz5);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx5, idy5, idz5);



		//----------Polarization 6------//
		mat_type = pt_glb->material_cell(idx6, idy6, idz6);
		//printf("mat_type6 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx6, idy6, idz6);

		if (mat_type == 0) {
			if ((idx6 < xS && true == if_PML_Xs) \
			|| (idx6 >= xE && true == if_PML_Xe) \
			|| (idy6 < yS && true == if_PML_Ys) \
			|| (idy6 >= yE && true == if_PML_Ye) \
			|| (idz6 < zS && true == if_PML_Zs) \
			|| (idz6 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx6, idy6, idz6); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx6, idy6, idz6); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx6, idy6, idz6); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx6, idy6, idz6);
		Jfy = Jfy + DJfy(idx6, idy6, idz6);
		Jfz = Jfz + DJfz(idx6, idy6, idz6);
		Jpx = Jpx + Jpx_n1_store(idx6, idy6, idz6) + Jpx_n2_store(idx6, idy6, idz6) + Jpx_n3_store(idx6, idy6, idz6);
		Jpy = Jpy + Jpy_n1_store(idx6, idy6, idz6) + Jpy_n2_store(idx6, idy6, idz6) + Jpy_n3_store(idx6, idy6, idz6);
		Jpz = Jpz + Jpz_n1_store(idx6, idy6, idz6) + Jpz_n2_store(idx6, idy6, idz6) + Jpz_n3_store(idx6, idy6, idz6);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx6, idy6, idz6);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx6, idy6, idz6);



		//----------Polarization 7------//
		mat_type = pt_glb->material_cell(idx7, idy7, idz7);
		//printf("mat_type7 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx7, idy7, idz7);

		if (mat_type == 0) {
			if ((idx7 < xS && true == if_PML_Xs) \
			|| (idx7 >= xE && true == if_PML_Xe) \
			|| (idy7 < yS && true == if_PML_Ys) \
			|| (idy7 >= yE && true == if_PML_Ye) \
			|| (idz7 < zS && true == if_PML_Zs) \
			|| (idz7 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx7, idy7, idz7); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx7, idy7, idz7); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx7, idy7, idz7); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx7, idy7, idz7);
		Jfy = Jfy + DJfy(idx7, idy7, idz7);
		Jfz = Jfz + DJfz(idx7, idy7, idz7);
		Jpx = Jpx + Jpx_n1_store(idx7, idy7, idz7) + Jpx_n2_store(idx7, idy7, idz7) + Jpx_n3_store(idx7, idy7, idz7);
		Jpy = Jpy + Jpy_n1_store(idx7, idy7, idz7) + Jpy_n2_store(idx7, idy7, idz7) + Jpy_n3_store(idx7, idy7, idz7);
		Jpz = Jpz + Jpz_n1_store(idx7, idy7, idz7) + Jpz_n2_store(idx7, idy7, idz7) + Jpz_n3_store(idx7, idy7, idz7);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx7, idy7, idz7);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx7, idy7, idz7);


		//----------Polarization 8------//
		mat_type = pt_glb->material_cell(idx8, idy8, idz8);
		//printf("mat_type8 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx8, idy8, idz8);
		
		if (mat_type == 0) {
			if ((idx8 < xS && true == if_PML_Xs) \
			|| (idx8 >= xE && true == if_PML_Xe) \
			|| (idy8 < yS && true == if_PML_Ys) \
			|| (idy8 >= yE && true == if_PML_Ye) \
			|| (idz8 < zS && true == if_PML_Zs) \
			|| (idz8 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(idx8, idy8, idz8); //RK
				dPy = dPy + pt_fe->dpy_glb_rk1(idx8, idy8, idz8); //RK
				dPz = dPz + pt_fe->dpz_glb_rk1(idx8, idy8, idz8); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx8, idy8, idz8);
		Jfy = Jfy + DJfy(idx8, idy8, idz8);
		Jfz = Jfz + DJfz(idx8, idy8, idz8);
		Jpx = Jpx + Jpx_n1_store(idx8, idy8, idz8) + Jpx_n2_store(idx8, idy8, idz8) + Jpx_n3_store(idx8, idy8, idz8);
		Jpy = Jpy + Jpy_n1_store(idx8, idy8, idz8) + Jpy_n2_store(idx8, idy8, idz8) + Jpy_n3_store(idx8, idy8, idz8);
		Jpz = Jpz + Jpz_n1_store(idx8, idy8, idz8) + Jpz_n2_store(idx8, idy8, idz8) + Jpz_n3_store(idx8, idy8, idz8);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx8, idy8, idz8);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx8, idy8, idz8);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jfx = Jfx / 8.;	Jfy = Jfy / 8.;	Jfz = Jfz / 8.;
		Jpx = Jpx / 8.;	Jpy = Jpy / 8.;	Jpz = Jpz / 8.;
		Jishex = Jishex / 8.;	Jishey = Jishey / 8.;
		dPx = dPx / P_count;	dPy = dPy / P_count;	dPz = dPz / P_count;

		cond11/=8.; cond12/=8.; cond13/=8.;
		cond21/=8.; cond22/=8.; cond23/=8.;
		cond31/=8.; cond32/=8.; cond33/=8.;

		er11/=8.; er12/=8.; er13/=8.;
		er21/=8.; er22/=8.; er23/=8.;
		er31/=8.; er32/=8.; er33/=8.;

		Denominator = -er13 * er22 * er31 + er12 * er23 * er31 + er13 * er21 * er32 - \
			er11 * er23 * er32 - er12 * er21 * er33 + er11 * er22 * er33;
		inverse_er11 = (-er23 * er32 + er22 * er33) / Denominator;
		inverse_er12 = (er13 * er32 - er12 * er33) / Denominator;
		inverse_er13 = (-er13 * er22 + er12 * er23) / Denominator;
		inverse_er21 = (er23 * er31 - er21 * er33) / Denominator;
		inverse_er22 = (-er13 * er31 + er11 * er33) / Denominator;
		inverse_er23 = (er13 * er21 - er11 * er23) / Denominator;
		inverse_er31 = (-er22 * er31 + er21 * er32) / Denominator;
		inverse_er32 = (er12 * er31 - er11 * er32) / Denominator;
		inverse_er33 = (-er12 * er21 + er11 * er22) / Denominator;

		isPML = false;

		if (((true == if_PML_Xs && i < xS) || (true == if_PML_Xe && i >= xE)) \
			|| ((true == if_PML_Ys && j < yS) || (true == if_PML_Ye && j >= yE)) \
			|| ((true == if_PML_Zs && k < zS) || (true == if_PML_Ze && k >= zE)))
		{
			isPML = true;

			if (i < nx) {
				DEx_em_t1(idx8, j, k) = pt_glb->dt / kappa_y_np1(0, j, 0) \
					* (dHzdy - dHydz - (sigma_y_np1(0, j, 0) / e0) * Dx_PML_store(idx8, j, k) \
						- Jfx - Jpx - Jishex - dPx / pt_glb->dt);
			}
			if (j < ny)
					DEy_em_t1(i, idy8, k) = pt_glb->dt / kappa_z_np1(0, 0, k) \
					* (dHxdz - dHzdx - (sigma_z_np1(0, 0, k) / e0) * Dy_PML_store(i, idy8, k) \
						- Jfy - Jpy - Jishey - dPy / pt_glb->dt);

			if (k < nz)
				DEz_em_t1(i, j, idz8) = pt_glb->dt / kappa_x_np1(i, 0, 0) \
					* (dHydx - dHxdy - (sigma_x_np1(i, 0, 0) / e0) * Dz_PML_store(i, j, idz8) \
					- Jfz - Jpz - dPz / pt_glb->dt);
		}

		if (true == isPML)
		{
			if (i < nx) {
				dDEx_em_rk1(idx8, j, k) = (pt_glb->dt / kappa_z_np1(0, 0, k) / e0) \
					* (inverse_er11*(kappa_x_n(idx8, 0, 0) * DEx_em_t1(idx8, j, k) / pt_glb->dt \
						+ sigma_x_n(idx8, 0, 0) * Dx_PML_store(idx8, j, k) / e0)\
						- (sigma_z_np1(0, 0, k) * DEx_em_store(idx8, j, k)));
			}

			if (j < ny) {
				dDEy_em_rk1(i, idy8, k) = (pt_glb->dt / kappa_x_np1(i, 0, 0) / e0) \
					* (inverse_er22*(kappa_y_n(0, idy8, 0) * DEy_em_t1(i, idy8, k) / pt_glb->dt \
						+ sigma_y_n(0, idy8, 0) * Dy_PML_store(i, idy8, k) / e0)\
						- (sigma_x_np1(i, 0, 0) * DEy_em_store(i, idy8, k)));
			}

			if (k < nz) {
				dDEz_em_rk1(i, j, idz8) = (pt_glb->dt / kappa_y_np1(0, j, 0) / e0) \
					* (inverse_er33*(kappa_z_n(0, 0, idz8) * DEz_em_t1(i, j, idz8) / pt_glb->dt \
						+ sigma_z_n(0, 0, idz8) * Dz_PML_store(i, j, idz8) / e0) \
						- (sigma_y_np1(0, j, 0) * DEz_em_store(i, j, idz8)));
			}
		}
		else
		{
			if (i < nx) {
				dDEx_em_rk1(idx8, j, k) = \
					pt_glb->dt / e0 * inverse_er11 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er12 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er13 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (j < ny) {

				dDEy_em_rk1(i, idy8, k) = \
					pt_glb->dt / e0 * inverse_er21 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er22 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er23 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (k < nz) {

				dDEz_em_rk1(i, j, idz8) = \
					pt_glb->dt / e0 * inverse_er31 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er32 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er33 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}
		}
	//-----------END X Y Z component-----------//
	}
}

void EMdynamic_system::get_dE_RK2() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8;
	long int idy1, idy2, idy3, idy4, idy5, idy6, idy7, idy8;
	long int idz1, idz2, idz3, idz4, idz5, idz6, idz7, idz8;

	double er11, er12, er13, er21, er22, er23, er31, er32, er33;

	double Denominator;
	double inverse_er11, inverse_er12, inverse_er13;
	double inverse_er21, inverse_er22, inverse_er23;
	double inverse_er31, inverse_er32, inverse_er33;


	double cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	bool isPML = false;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_half(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK2();
	}

#pragma acc wait




	//-----------X Y Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1, idx2, idy2, idz2, idx3, idy3, idz3, idx4, idy4, idz4,\
idx5, idy5, idz5, idx6, idy6, idz6, idx7, idy7, idz7, idx8, idy8, idz8,\
er11, er12, er13, er21, er22, er23, er31, er32, er33,\
cond11, con12, cond13, cond21, con22, cond23, cond31, con32, cond33,\
Denominator,\
inverse_er11, inverse_er12, inverse_er13,\
inverse_er21, inverse_er22, inverse_er23,\
inverse_er31, inverse_er32, inverse_er33,\
cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,isPML,  \
P_count) default(present)async(1)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);


		if ((i == 0 || i == nx) && periodicX_EM == false) {
			if (if_PML_Xe == false && if_PML_Xs == false) {
				continue;
			}
			else {
				idx1 = i; idx2 = i; idx3 = i; idx4 = i; idx5 = i; idx6 = i; idx7 = i; idx8 = i;	
			}
		}
		else if ((i == 0 || i == nx) && periodicX_EM == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = nx - 1; idx4 = nx - 1; idx5 = 0; idx6 = 0; idx7 = 0; idx8 = 0;
		}
		else if (i < nx + 1) {
			idx1 = i - 1; idx2 = i - 1; idx3 = i - 1; idx4 = i - 1; idx5 = i; idx6 = i; idx7 = i; idx8 = i;
		}
		if ((j == 0 || j == ny) && periodicY_EM == false) {
			if (if_PML_Ye == false && if_PML_Ys == false) {
				continue;
			}
			else {
			idy1 = j; idy2 = j; idy3 = j; idy4 = j; idy5 = j; idy6 = j; idy7 = j; idy8 = j;
			}
		}
		else if ((j == 0 || j == ny) && periodicY_EM == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0; idy5 = ny - 1; idy6 = ny - 1; idy7 = 0; idy8 = 0;
		}
		else if (j < ny + 1) {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j; idy5 = j - 1; idy6 = j - 1; idy7 = j; idy8 = j;
		}
		if ((k == 0 || k == nz) && periodicZ_EM == false) {
			if (if_PML_Ze == false && if_PML_Zs == false) {
				continue;
			}
			else {
			idz1 = k; idz2 = k; idz3 = k; idz4 = k; idz5 = k; idz6 = k; idz7 = k; idz8 = k;
			}
		}
		else if ((k == 0 || k == nz) && periodicZ_EM == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0; idz5 = nz - 1; idz6 = 0; idz7 = nz - 1; idz8 = 0;
		}
		else if (k < nz + 1) {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k; idz5 = k - 1; idz6 = k; idz7 = k - 1; idz8 = k;
		}

		dHzdy = (DHz_em_store(idx5, idy3, idz4) - DHz_em_store(idx5, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx5, idy3, idz1)) / dz;
		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.;

		dHxdz = (DHx_em_store(idx5, idy7, idz2) - DHx_em_store(idx5, idy7, idz1)) / dz;
		dHzdx = (DHz_em_store(idx5, idy7, idz2) - DHz_em_store(idx1, idy7, idz2)) / dx;
		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.;

		dHydx = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx1, idy3, idz2)) / dx;
		dHxdy = (DHx_em_store(idx5, idy3, idz2) - DHx_em_store(idx5, idy1, idz2)) / dy;
		Jfz = 0.; Jpz = 0.;
		dPz = 0.;

		cond11 = 0.; cond12 = 0.; cond13 = 0.; cond21 = 0.; cond22 = 0.; cond23 = 0.; cond31 = 0.; cond32 = 0.; cond33 = 0.;
		er11 = 0.; er12 = 0.; er13 = 0.; er21 = 0.; er22 = 0.; er23 = 0.; er31 = 0.; er32 = 0.; er33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;



		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			if ((idx1 < xS && true == if_PML_Xs) \
			|| (idx1 >= xE && true == if_PML_Xe) \
			|| (idy1 < yS && true == if_PML_Ys) \
			|| (idy1 >= yE && true == if_PML_Ye) \
			|| (idz1 < zS && true == if_PML_Zs) \
			|| (idz1 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx1, idy1, idz1); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx1, idy1, idz1); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx1, idy1, idz1);
		Jfy = Jfy + DJfy(idx1, idy1, idz1);
		Jfz = Jfz + DJfz(idx1, idy1, idz1);
		Jpx = Jpx + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jpy = Jpy + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jpz = Jpz + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx1, idy1, idz1);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx1, idy1, idz1);


		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			if ((idx2 < xS && true == if_PML_Xs) \
				|| (idx2 >= xE && true == if_PML_Xe) \
				|| (idy2 < yS && true == if_PML_Ys) \
				|| (idy2 >= yE && true == if_PML_Ye) \
				|| (idz2 < zS && true == if_PML_Zs) \
				|| (idz2 >= zE && true == if_PML_Ze))
			{
				isPML = true;
			}
			else
				isPML = false;

			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx2, idy2, idz2); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx2, idy2, idz2); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx2, idy2, idz2);
		Jfy = Jfy + DJfy(idx2, idy2, idz2);
		Jfz = Jfz + DJfz(idx2, idy2, idz2);
		Jpx = Jpx + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jpy = Jpy + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jpz = Jpz + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx2, idy2, idz2);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx2, idy2, idz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		////printf("mat_type3 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx3, idy3, idz3);
		if (mat_type == 0) {
			if ((idx3 < xS && true == if_PML_Xs) \
			|| (idx3 >= xE && true == if_PML_Xe) \
			|| (idy3 < yS && true == if_PML_Ys) \
			|| (idy3 >= yE && true == if_PML_Ye) \
			|| (idz3 < zS && true == if_PML_Zs) \
			|| (idz3 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx3, idy3, idz3); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx3, idy3, idz3); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx3, idy3, idz3);
		Jfy = Jfy + DJfy(idx3, idy3, idz3);
		Jfz = Jfz + DJfz(idx3, idy3, idz3);
		Jpx = Jpx + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jpy = Jpy + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jpz = Jpz + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx3, idy3, idz3);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx3, idy3, idz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		//printf("mat_type4 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx4, idy4, idz4);
		if (mat_type == 0) {
			if ((idx4 < xS && true == if_PML_Xs) \
			|| (idx4 >= xE && true == if_PML_Xe) \
			|| (idy4 < yS && true == if_PML_Ys) \
			|| (idy4 >= yE && true == if_PML_Ye) \
			|| (idz4 < zS && true == if_PML_Zs) \
			|| (idz4 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx4, idy4, idz4); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx4, idy4, idz4); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx4, idy4, idz4);
		Jfy = Jfy + DJfy(idx4, idy4, idz4);
		Jfz = Jfz + DJfz(idx4, idy4, idz4);
		Jpx = Jpx + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jpy = Jpy + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jpz = Jpz + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx4, idy4, idz4);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx4, idy4, idz4);


		//----------Polarization 5------//
		mat_type = pt_glb->material_cell(idx5, idy5, idz5);
		if (mat_type == 0) {
			if ((idx5 < xS && true == if_PML_Xs) \
			|| (idx5 >= xE && true == if_PML_Xe) \
			|| (idy5 < yS && true == if_PML_Ys) \
			|| (idy5 >= yE && true == if_PML_Ye) \
			|| (idz5 < zS && true == if_PML_Zs) \
			|| (idz5 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx5, idy5, idz5); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx5, idy5, idz5); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx5, idy5, idz5); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx5, idy5, idz5);
		Jfy = Jfy + DJfy(idx5, idy5, idz5);
		Jfz = Jfz + DJfz(idx5, idy5, idz5);
		Jpx = Jpx + Jpx_n1_store(idx5, idy5, idz5) + Jpx_n2_store(idx5, idy5, idz5) + Jpx_n3_store(idx5, idy5, idz5);
		Jpy = Jpy + Jpy_n1_store(idx5, idy5, idz5) + Jpy_n2_store(idx5, idy5, idz5) + Jpy_n3_store(idx5, idy5, idz5);
		Jpz = Jpz + Jpz_n1_store(idx5, idy5, idz5) + Jpz_n2_store(idx5, idy5, idz5) + Jpz_n3_store(idx5, idy5, idz5);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx5, idy5, idz5);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx5, idy5, idz5);



		//----------Polarization 6------//
		mat_type = pt_glb->material_cell(idx6, idy6, idz6);
		//printf("mat_type6 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx6, idy6, idz6);

		if (mat_type == 0) {
			if ((idx6 < xS && true == if_PML_Xs) \
			|| (idx6 >= xE && true == if_PML_Xe) \
			|| (idy6 < yS && true == if_PML_Ys) \
			|| (idy6 >= yE && true == if_PML_Ye) \
			|| (idz6 < zS && true == if_PML_Zs) \
			|| (idz6 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx6, idy6, idz6); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx6, idy6, idz6); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx6, idy6, idz6); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx6, idy6, idz6);
		Jfy = Jfy + DJfy(idx6, idy6, idz6);
		Jfz = Jfz + DJfz(idx6, idy6, idz6);
		Jpx = Jpx + Jpx_n1_store(idx6, idy6, idz6) + Jpx_n2_store(idx6, idy6, idz6) + Jpx_n3_store(idx6, idy6, idz6);
		Jpy = Jpy + Jpy_n1_store(idx6, idy6, idz6) + Jpy_n2_store(idx6, idy6, idz6) + Jpy_n3_store(idx6, idy6, idz6);
		Jpz = Jpz + Jpz_n1_store(idx6, idy6, idz6) + Jpz_n2_store(idx6, idy6, idz6) + Jpz_n3_store(idx6, idy6, idz6);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx6, idy6, idz6);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx6, idy6, idz6);



		//----------Polarization 7------//
		mat_type = pt_glb->material_cell(idx7, idy7, idz7);
		//printf("mat_type7 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx7, idy7, idz7);

		if (mat_type == 0) {
			if ((idx7 < xS && true == if_PML_Xs) \
			|| (idx7 >= xE && true == if_PML_Xe) \
			|| (idy7 < yS && true == if_PML_Ys) \
			|| (idy7 >= yE && true == if_PML_Ye) \
			|| (idz7 < zS && true == if_PML_Zs) \
			|| (idz7 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx7, idy7, idz7); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx7, idy7, idz7); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx7, idy7, idz7); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx7, idy7, idz7);
		Jfy = Jfy + DJfy(idx7, idy7, idz7);
		Jfz = Jfz + DJfz(idx7, idy7, idz7);
		Jpx = Jpx + Jpx_n1_store(idx7, idy7, idz7) + Jpx_n2_store(idx7, idy7, idz7) + Jpx_n3_store(idx7, idy7, idz7);
		Jpy = Jpy + Jpy_n1_store(idx7, idy7, idz7) + Jpy_n2_store(idx7, idy7, idz7) + Jpy_n3_store(idx7, idy7, idz7);
		Jpz = Jpz + Jpz_n1_store(idx7, idy7, idz7) + Jpz_n2_store(idx7, idy7, idz7) + Jpz_n3_store(idx7, idy7, idz7);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx7, idy7, idz7);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx7, idy7, idz7);


		//----------Polarization 8------//
		mat_type = pt_glb->material_cell(idx8, idy8, idz8);
		
		if (mat_type == 0) {
			if ((idx8 < xS && true == if_PML_Xs) \
			|| (idx8 >= xE && true == if_PML_Xe) \
			|| (idy8 < yS && true == if_PML_Ys) \
			|| (idy8 >= yE && true == if_PML_Ye) \
			|| (idz8 < zS && true == if_PML_Zs) \
			|| (idz8 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(idx8, idy8, idz8); //RK
				dPy = dPy + pt_fe->dpy_glb_rk2(idx8, idy8, idz8); //RK
				dPz = dPz + pt_fe->dpz_glb_rk2(idx8, idy8, idz8); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx8, idy8, idz8);
		Jfy = Jfy + DJfy(idx8, idy8, idz8);
		Jfz = Jfz + DJfz(idx8, idy8, idz8);
		Jpx = Jpx + Jpx_n1_store(idx8, idy8, idz8) + Jpx_n2_store(idx8, idy8, idz8) + Jpx_n3_store(idx8, idy8, idz8);
		Jpy = Jpy + Jpy_n1_store(idx8, idy8, idz8) + Jpy_n2_store(idx8, idy8, idz8) + Jpy_n3_store(idx8, idy8, idz8);
		Jpz = Jpz + Jpz_n1_store(idx8, idy8, idz8) + Jpz_n2_store(idx8, idy8, idz8) + Jpz_n3_store(idx8, idy8, idz8);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx8, idy8, idz8);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx8, idy8, idz8);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jfx = Jfx / 8.;	Jfy = Jfy / 8.;	Jfz = Jfz / 8.;
		Jpx = Jpx / 8.;	Jpy = Jpy / 8.;	Jpz = Jpz / 8.;
		Jishex = Jishex / 8.;	Jishey = Jishey / 8.;
		dPx = dPx / P_count;	dPy = dPy / P_count;	dPz = dPz / P_count;

		cond11 = cond11 / 8.; cond12 = cond12 / 8.; cond13 = cond13 / 8.;
		cond21 = cond21 / 8.; cond22 = cond22 / 8.; cond23 = cond23 / 8.;
		cond31 = cond31 / 8.; cond32 = cond32 / 8.; cond33 = cond33 / 8.;

		er11 = er11 / 8.; er12 = er12 / 8.; er13 = er13 / 8.;
		er21 = er21 / 8.; er22 = er22 / 8.; er23 = er23 / 8.;
		er31 = er31 / 8.; er32 = er32 / 8.; er33 = er33 / 8.;

		//printf("%ld\t%ld\t%ld\t%le\t%le\t%le\n", i, j, k, er11, er22, er33);


		Denominator = -er13 * er22 * er31 + er12 * er23 * er31 + er13 * er21 * er32 - \
			er11 * er23 * er32 - er12 * er21 * er33 + er11 * er22 * er33;
		inverse_er11 = (-er23 * er32 + er22 * er33) / Denominator;
		inverse_er12 = (er13 * er32 - er12 * er33) / Denominator;
		inverse_er13 = (-er13 * er22 + er12 * er23) / Denominator;
		inverse_er21 = (er23 * er31 - er21 * er33) / Denominator;
		inverse_er22 = (-er13 * er31 + er11 * er33) / Denominator;
		inverse_er23 = (er13 * er21 - er11 * er23) / Denominator;
		inverse_er31 = (-er22 * er31 + er21 * er32) / Denominator;
		inverse_er32 = (er12 * er31 - er11 * er32) / Denominator;
		inverse_er33 = (-er12 * er21 + er11 * er22) / Denominator;

		isPML = false;

		if (((true == if_PML_Xs && i < xS) || (true == if_PML_Xe && i >= xE)) \
			|| ((true == if_PML_Ys && j < yS) || (true == if_PML_Ye && j >= yE)) \
			|| ((true == if_PML_Zs && k < zS) || (true == if_PML_Ze && k >= zE)))
		{
			isPML = true;

			if (i < nx)
				DEx_em_t2(idx8, j, k) = pt_glb->dt / kappa_y_np1(0, j, 0) \
					* (dHzdy - dHydz - (sigma_y_np1(0, j, 0) / e0) * Dx_PML_store(idx8, j, k) \
					- Jfx - Jpx - Jishex - dPx / pt_glb->dt);

			if (j < ny)
				DEy_em_t2(i, idy8, k) = pt_glb->dt / kappa_z_np1(0, 0, k) \
					* (dHxdz - dHzdx - (sigma_z_np1(0, 0, k) / e0) * Dy_PML_store(i, idy8, k) \
						- Jfy - Jpy - Jishey - dPy / pt_glb->dt);

			if (k < nz)
				DEz_em_t2(i, j, idz8) = pt_glb->dt / kappa_x_np1(i, 0, 0) \
					* (dHydx - dHxdy - (sigma_x_np1(i, 0, 0) / e0) * Dz_PML_store(i, j, idz8) \
					- Jfz - Jpz - dPz / pt_glb->dt);
		}

		if (true == isPML)
		{
			if (i < nx) {
				dDEx_em_rk2(idx8, j, k) = (pt_glb->dt / kappa_z_np1(0, 0, k) / e0) \
					* (inverse_er11*(kappa_x_n(idx8, 0, 0) * DEx_em_t2(idx8, j, k) / pt_glb->dt \
						+ sigma_x_n(idx8, 0, 0) * Dx_PML_store(idx8, j, k) / e0)\
						- (sigma_z_np1(0, 0, k) * DEx_em_store(idx8, j, k)));
			}

			if (j < ny) {
				dDEy_em_rk2(i, idy8, k) = (pt_glb->dt / kappa_x_np1(i, 0, 0) / e0) \
					* (inverse_er22*(kappa_y_n(0, idy8, 0) * DEy_em_t2(i, idy8, k) / pt_glb->dt \
						+ sigma_y_n(0, idy8, 0) * Dy_PML_store(i, idy8, k) / e0)\
						- (sigma_x_np1(i, 0, 0) * DEy_em_store(i, idy8, k)));
			}

			if (k < nz) {
				dDEz_em_rk2(i, j, idz8) = (pt_glb->dt / kappa_y_np1(0, j, 0) / e0) \
					* (inverse_er33*(kappa_z_n(0, 0, idz8) * DEz_em_t2(i, j, idz8) / pt_glb->dt \
						+ sigma_z_n(0, 0, idz8) * Dz_PML_store(i, j, idz8) / e0) \
						- (sigma_y_np1(0, j, 0) * DEz_em_store(i, j, idz8)));
			}
		}
		else
		{
			if (i < nx) {
				dDEx_em_rk2(idx8, j, k) = \
					pt_glb->dt / e0 * inverse_er11 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er12 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er13 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (j < ny) {
				dDEy_em_rk2(i, idy8, k) = \
					pt_glb->dt / e0 * inverse_er21 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er22 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er23 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (k < nz) {
				dDEz_em_rk2(i, j, idz8) = \
					pt_glb->dt / e0 * inverse_er31 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er32 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er33 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}
		}
	}
	//-----------END X Y Z component-----------//
}


void EMdynamic_system::get_dE_RK3() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8;
	long int idy1, idy2, idy3, idy4, idy5, idy6, idy7, idy8;
	long int idz1, idz2, idz3, idz4, idz5, idz6, idz7, idz8;

	double er11, er12, er13, er21, er22, er23, er31, er32, er33;

	double Denominator;
	double inverse_er11, inverse_er12, inverse_er13;
	double inverse_er21, inverse_er22, inverse_er23;
	double inverse_er31, inverse_er32, inverse_er33;


	double cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	bool isPML = false;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_half(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK3();
	}

#pragma acc wait

	//-----------X Y Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1, idx2, idy2, idz2, idx3, idy3, idz3, idx4, idy4, idz4,\
idx5, idy5, idz5, idx6, idy6, idz6, idx7, idy7, idz7, idx8, idy8, idz8,\
er11, er12, er13, er21, er22, er23, er31, er32, er33,\
cond11, con12, cond13, cond21, con22, cond23, cond31, con32, cond33,\
Denominator,\
inverse_er11, inverse_er12, inverse_er13,\
inverse_er21, inverse_er22, inverse_er23,\
inverse_er31, inverse_er32, inverse_er33,\
cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,isPML,  \
P_count) default(present)async(1)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);


		if ((i == 0 || i == nx) && periodicX_EM == false) {
			if (if_PML_Xe == false && if_PML_Xs == false) {
				continue;
			}
			else {
				idx1 = i; idx2 = i; idx3 = i; idx4 = i; idx5 = i; idx6 = i; idx7 = i; idx8 = i;	
			}
		}
		else if ((i == 0 || i == nx) && periodicX_EM == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = nx - 1; idx4 = nx - 1; idx5 = 0; idx6 = 0; idx7 = 0; idx8 = 0;
		}
		else if (i < nx + 1) {
			idx1 = i - 1; idx2 = i - 1; idx3 = i - 1; idx4 = i - 1; idx5 = i; idx6 = i; idx7 = i; idx8 = i;
		}
		if ((j == 0 || j == ny) && periodicY_EM == false) {
			if (if_PML_Ye == false && if_PML_Ys == false) {
				continue;
			}
			else {
			idy1 = j; idy2 = j; idy3 = j; idy4 = j; idy5 = j; idy6 = j; idy7 = j; idy8 = j;
			}
		}
		else if ((j == 0 || j == ny) && periodicY_EM == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0; idy5 = ny - 1; idy6 = ny - 1; idy7 = 0; idy8 = 0;
		}
		else if (j < ny + 1) {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j; idy5 = j - 1; idy6 = j - 1; idy7 = j; idy8 = j;
		}
		if ((k == 0 || k == nz) && periodicZ_EM == false) {
			if (if_PML_Ze == false && if_PML_Zs == false) {
				continue;
			}
			else {
			idz1 = k; idz2 = k; idz3 = k; idz4 = k; idz5 = k; idz6 = k; idz7 = k; idz8 = k;
			}
		}
		else if ((k == 0 || k == nz) && periodicZ_EM == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0; idz5 = nz - 1; idz6 = 0; idz7 = nz - 1; idz8 = 0;
		}
		else if (k < nz + 1) {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k; idz5 = k - 1; idz6 = k; idz7 = k - 1; idz8 = k;
		}

		dHzdy = (DHz_em_store(idx5, idy3, idz4) - DHz_em_store(idx5, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx5, idy3, idz1)) / dz;
		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.;

		dHxdz = (DHx_em_store(idx5, idy7, idz2) - DHx_em_store(idx5, idy7, idz1)) / dz;
		dHzdx = (DHz_em_store(idx5, idy7, idz2) - DHz_em_store(idx1, idy7, idz2)) / dx;
		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.;

		dHydx = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx1, idy3, idz2)) / dx;
		dHxdy = (DHx_em_store(idx5, idy3, idz2) - DHx_em_store(idx5, idy1, idz2)) / dy;
		Jfz = 0.; Jpz = 0.;
		dPz = 0.;

		cond11 = 0.; cond12 = 0.; cond13 = 0.; cond21 = 0.; cond22 = 0.; cond23 = 0.; cond31 = 0.; cond32 = 0.; cond33 = 0.;
		er11 = 0.; er12 = 0.; er13 = 0.; er21 = 0.; er22 = 0.; er23 = 0.; er31 = 0.; er32 = 0.; er33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;



		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			if ((idx1 < xS && true == if_PML_Xs) \
			|| (idx1 >= xE && true == if_PML_Xe) \
			|| (idy1 < yS && true == if_PML_Ys) \
			|| (idy1 >= yE && true == if_PML_Ye) \
			|| (idz1 < zS && true == if_PML_Zs) \
			|| (idz1 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx1, idy1, idz1); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx1, idy1, idz1); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx1, idy1, idz1);
		Jfy = Jfy + DJfy(idx1, idy1, idz1);
		Jfz = Jfz + DJfz(idx1, idy1, idz1);
		Jpx = Jpx + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jpy = Jpy + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jpz = Jpz + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx1, idy1, idz1);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx1, idy1, idz1);


		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			if ((idx2 < xS && true == if_PML_Xs) \
				|| (idx2 >= xE && true == if_PML_Xe) \
				|| (idy2 < yS && true == if_PML_Ys) \
				|| (idy2 >= yE && true == if_PML_Ye) \
				|| (idz2 < zS && true == if_PML_Zs) \
				|| (idz2 >= zE && true == if_PML_Ze))
			{
				isPML = true;
			}
			else
				isPML = false;

			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx2, idy2, idz2); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx2, idy2, idz2); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx2, idy2, idz2);
		Jfy = Jfy + DJfy(idx2, idy2, idz2);
		Jfz = Jfz + DJfz(idx2, idy2, idz2);
		Jpx = Jpx + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jpy = Jpy + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jpz = Jpz + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx2, idy2, idz2);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx2, idy2, idz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		////printf("mat_type3 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx3, idy3, idz3);
		if (mat_type == 0) {
			if ((idx3 < xS && true == if_PML_Xs) \
			|| (idx3 >= xE && true == if_PML_Xe) \
			|| (idy3 < yS && true == if_PML_Ys) \
			|| (idy3 >= yE && true == if_PML_Ye) \
			|| (idz3 < zS && true == if_PML_Zs) \
			|| (idz3 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx3, idy3, idz3); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx3, idy3, idz3); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx3, idy3, idz3);
		Jfy = Jfy + DJfy(idx3, idy3, idz3);
		Jfz = Jfz + DJfz(idx3, idy3, idz3);
		Jpx = Jpx + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jpy = Jpy + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jpz = Jpz + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx3, idy3, idz3);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx3, idy3, idz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		//printf("mat_type4 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx4, idy4, idz4);
		if (mat_type == 0) {
			if ((idx4 < xS && true == if_PML_Xs) \
			|| (idx4 >= xE && true == if_PML_Xe) \
			|| (idy4 < yS && true == if_PML_Ys) \
			|| (idy4 >= yE && true == if_PML_Ye) \
			|| (idz4 < zS && true == if_PML_Zs) \
			|| (idz4 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx4, idy4, idz4); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx4, idy4, idz4); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx4, idy4, idz4);
		Jfy = Jfy + DJfy(idx4, idy4, idz4);
		Jfz = Jfz + DJfz(idx4, idy4, idz4);
		Jpx = Jpx + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jpy = Jpy + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jpz = Jpz + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx4, idy4, idz4);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx4, idy4, idz4);


		//----------Polarization 5------//
		mat_type = pt_glb->material_cell(idx5, idy5, idz5);
		if (mat_type == 0) {
			if ((idx5 < xS && true == if_PML_Xs) \
			|| (idx5 >= xE && true == if_PML_Xe) \
			|| (idy5 < yS && true == if_PML_Ys) \
			|| (idy5 >= yE && true == if_PML_Ye) \
			|| (idz5 < zS && true == if_PML_Zs) \
			|| (idz5 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx5, idy5, idz5); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx5, idy5, idz5); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx5, idy5, idz5); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx5, idy5, idz5);
		Jfy = Jfy + DJfy(idx5, idy5, idz5);
		Jfz = Jfz + DJfz(idx5, idy5, idz5);
		Jpx = Jpx + Jpx_n1_store(idx5, idy5, idz5) + Jpx_n2_store(idx5, idy5, idz5) + Jpx_n3_store(idx5, idy5, idz5);
		Jpy = Jpy + Jpy_n1_store(idx5, idy5, idz5) + Jpy_n2_store(idx5, idy5, idz5) + Jpy_n3_store(idx5, idy5, idz5);
		Jpz = Jpz + Jpz_n1_store(idx5, idy5, idz5) + Jpz_n2_store(idx5, idy5, idz5) + Jpz_n3_store(idx5, idy5, idz5);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx5, idy5, idz5);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx5, idy5, idz5);



		//----------Polarization 6------//
		mat_type = pt_glb->material_cell(idx6, idy6, idz6);
		//printf("mat_type6 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx6, idy6, idz6);

		if (mat_type == 0) {
			if ((idx6 < xS && true == if_PML_Xs) \
			|| (idx6 >= xE && true == if_PML_Xe) \
			|| (idy6 < yS && true == if_PML_Ys) \
			|| (idy6 >= yE && true == if_PML_Ye) \
			|| (idz6 < zS && true == if_PML_Zs) \
			|| (idz6 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx6, idy6, idz6); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx6, idy6, idz6); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx6, idy6, idz6); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx6, idy6, idz6);
		Jfy = Jfy + DJfy(idx6, idy6, idz6);
		Jfz = Jfz + DJfz(idx6, idy6, idz6);
		Jpx = Jpx + Jpx_n1_store(idx6, idy6, idz6) + Jpx_n2_store(idx6, idy6, idz6) + Jpx_n3_store(idx6, idy6, idz6);
		Jpy = Jpy + Jpy_n1_store(idx6, idy6, idz6) + Jpy_n2_store(idx6, idy6, idz6) + Jpy_n3_store(idx6, idy6, idz6);
		Jpz = Jpz + Jpz_n1_store(idx6, idy6, idz6) + Jpz_n2_store(idx6, idy6, idz6) + Jpz_n3_store(idx6, idy6, idz6);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx6, idy6, idz6);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx6, idy6, idz6);



		//----------Polarization 7------//
		mat_type = pt_glb->material_cell(idx7, idy7, idz7);
		//printf("mat_type7 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx7, idy7, idz7);

		if (mat_type == 0) {
			if ((idx7 < xS && true == if_PML_Xs) \
			|| (idx7 >= xE && true == if_PML_Xe) \
			|| (idy7 < yS && true == if_PML_Ys) \
			|| (idy7 >= yE && true == if_PML_Ye) \
			|| (idz7 < zS && true == if_PML_Zs) \
			|| (idz7 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx7, idy7, idz7); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx7, idy7, idz7); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx7, idy7, idz7); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx7, idy7, idz7);
		Jfy = Jfy + DJfy(idx7, idy7, idz7);
		Jfz = Jfz + DJfz(idx7, idy7, idz7);
		Jpx = Jpx + Jpx_n1_store(idx7, idy7, idz7) + Jpx_n2_store(idx7, idy7, idz7) + Jpx_n3_store(idx7, idy7, idz7);
		Jpy = Jpy + Jpy_n1_store(idx7, idy7, idz7) + Jpy_n2_store(idx7, idy7, idz7) + Jpy_n3_store(idx7, idy7, idz7);
		Jpz = Jpz + Jpz_n1_store(idx7, idy7, idz7) + Jpz_n2_store(idx7, idy7, idz7) + Jpz_n3_store(idx7, idy7, idz7);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx7, idy7, idz7);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx7, idy7, idz7);


		//----------Polarization 8------//
		mat_type = pt_glb->material_cell(idx8, idy8, idz8);
		//printf("mat_type8 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx8, idy8, idz8);
		
		if (mat_type == 0) {
			if ((idx8 < xS && true == if_PML_Xs) \
			|| (idx8 >= xE && true == if_PML_Xe) \
			|| (idy8 < yS && true == if_PML_Ys) \
			|| (idy8 >= yE && true == if_PML_Ye) \
			|| (idz8 < zS && true == if_PML_Zs) \
			|| (idz8 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(idx8, idy8, idz8); //RK
				dPy = dPy + pt_fe->dpy_glb_rk3(idx8, idy8, idz8); //RK
				dPz = dPz + pt_fe->dpz_glb_rk3(idx8, idy8, idz8); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx8, idy8, idz8);
		Jfy = Jfy + DJfy(idx8, idy8, idz8);
		Jfz = Jfz + DJfz(idx8, idy8, idz8);
		Jpx = Jpx + Jpx_n1_store(idx8, idy8, idz8) + Jpx_n2_store(idx8, idy8, idz8) + Jpx_n3_store(idx8, idy8, idz8);
		Jpy = Jpy + Jpy_n1_store(idx8, idy8, idz8) + Jpy_n2_store(idx8, idy8, idz8) + Jpy_n3_store(idx8, idy8, idz8);
		Jpz = Jpz + Jpz_n1_store(idx8, idy8, idz8) + Jpz_n2_store(idx8, idy8, idz8) + Jpz_n3_store(idx8, idy8, idz8);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx8, idy8, idz8);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx8, idy8, idz8);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jfx = Jfx / 8.;	Jfy = Jfy / 8.;	Jfz = Jfz / 8.;
		Jpx = Jpx / 8.;	Jpy = Jpy / 8.;	Jpz = Jpz / 8.;
		Jishex = Jishex / 8.;	Jishey = Jishey / 8.;
		dPx = dPx / P_count;	dPy = dPy / P_count;	dPz = dPz / P_count;

		cond11 = cond11 / 8.; cond12 = cond12 / 8.; cond13 = cond13 / 8.;
		cond21 = cond21 / 8.; cond22 = cond22 / 8.; cond23 = cond23 / 8.;
		cond31 = cond31 / 8.; cond32 = cond32 / 8.; cond33 = cond33 / 8.;

		er11 = er11 / 8.; er12 = er12 / 8.; er13 = er13 / 8.;
		er21 = er21 / 8.; er22 = er22 / 8.; er23 = er23 / 8.;
		er31 = er31 / 8.; er32 = er32 / 8.; er33 = er33 / 8.;


		Denominator = -er13 * er22 * er31 + er12 * er23 * er31 + er13 * er21 * er32 - \
			er11 * er23 * er32 - er12 * er21 * er33 + er11 * er22 * er33;
		inverse_er11 = (-er23 * er32 + er22 * er33) / Denominator;
		inverse_er12 = (er13 * er32 - er12 * er33) / Denominator;
		inverse_er13 = (-er13 * er22 + er12 * er23) / Denominator;
		inverse_er21 = (er23 * er31 - er21 * er33) / Denominator;
		inverse_er22 = (-er13 * er31 + er11 * er33) / Denominator;
		inverse_er23 = (er13 * er21 - er11 * er23) / Denominator;
		inverse_er31 = (-er22 * er31 + er21 * er32) / Denominator;
		inverse_er32 = (er12 * er31 - er11 * er32) / Denominator;
		inverse_er33 = (-er12 * er21 + er11 * er22) / Denominator;

		isPML = false;

		if (((true == if_PML_Xs && i < xS) || (true == if_PML_Xe && i >= xE)) \
			|| ((true == if_PML_Ys && j < yS) || (true == if_PML_Ye && j >= yE)) \
			|| ((true == if_PML_Zs && k < zS) || (true == if_PML_Ze && k >= zE)))
		{
			isPML = true;

			if (i < nx)
				DEx_em_t3(idx8, j, k) = pt_glb->dt / kappa_y_np1(0, j, 0) \
					* (dHzdy - dHydz - (sigma_y_np1(0, j, 0) / e0) * Dx_PML_store(idx8, j, k) \
					- Jfx - Jpx - Jishex - dPx / pt_glb->dt);

			if (j < ny)
				DEy_em_t3(i, idy8, k) = pt_glb->dt / kappa_z_np1(0, 0, k) \
					* (dHxdz - dHzdx - (sigma_z_np1(0, 0, k) / e0) * Dy_PML_store(i, idy8, k) \
						- Jfy - Jpy - Jishey - dPy / pt_glb->dt);

			if (k < nz)
				DEz_em_t3(i, j, idz8) = pt_glb->dt / kappa_x_np1(i, 0, 0) \
					* (dHydx - dHxdy - (sigma_x_np1(i, 0, 0) / e0) * Dz_PML_store(i, j, idz8) \
					- Jfz - Jpz - dPz / pt_glb->dt);
		}

		if (true == isPML)
		{
			if (i < nx) {
				dDEx_em_rk3(idx8, j, k) = (pt_glb->dt / kappa_z_np1(0, 0, k) / e0) \
					* (inverse_er11*(kappa_x_n(idx8, 0, 0) * DEx_em_t3(idx8, j, k) / pt_glb->dt \
						+ sigma_x_n(idx8, 0, 0) * Dx_PML_store(idx8, j, k) / e0)\
						- (sigma_z_np1(0, 0, k) * DEx_em_store(idx8, j, k)));
			}

			if (j < ny) {
				dDEy_em_rk3(i, idy8, k) = (pt_glb->dt / kappa_x_np1(i, 0, 0) / e0) \
					* (inverse_er22*(kappa_y_n(0, idy8, 0) * DEy_em_t3(i, idy8, k) / pt_glb->dt \
						+ sigma_y_n(0, idy8, 0) * Dy_PML_store(i, idy8, k) / e0)\
						- (sigma_x_np1(i, 0, 0) * DEy_em_store(i, idy8, k)));
			}

			if (k < nz) {
				dDEz_em_rk3(i, j, idz8) = (pt_glb->dt / kappa_y_np1(0, j, 0) / e0) \
					* (inverse_er33*(kappa_z_n(0, 0, idz8) * DEz_em_t3(i, j, idz8) / pt_glb->dt \
						+ sigma_z_n(0, 0, idz8) * Dz_PML_store(i, j, idz8) / e0) \
						- (sigma_y_np1(0, j, 0) * DEz_em_store(i, j, idz8)));
			}
		}
		else
		{
			if (i < nx) {
				dDEx_em_rk3(idx8, j, k) = \
					pt_glb->dt / e0 * inverse_er11 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er12 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er13 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (j < ny) {
				dDEy_em_rk3(i, idy8, k) = \
					pt_glb->dt / e0 * inverse_er21 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er22 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er23 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (k < nz) {
				dDEz_em_rk3(i, j, idz8) = \
					pt_glb->dt / e0 * inverse_er31 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er32 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er33 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}
		}
	}
	//-----------END X Y Z component-----------//

}

void EMdynamic_system::get_dE_RK4() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8;
	long int idy1, idy2, idy3, idy4, idy5, idy6, idy7, idy8;
	long int idz1, idz2, idz3, idz4, idz5, idz6, idz7, idz8;

	double er11, er12, er13, er21, er22, er23, er31, er32, er33;

	double Denominator;
	double inverse_er11, inverse_er12, inverse_er13;
	double inverse_er21, inverse_er22, inverse_er23;
	double inverse_er31, inverse_er32, inverse_er33;


	double cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	bool isPML = false;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_full(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK4();
	}

#pragma acc wait



	//-----------X Y Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1, idx2, idy2, idz2, idx3, idy3, idz3, idx4, idy4, idz4,\
idx5, idy5, idz5, idx6, idy6, idz6, idx7, idy7, idz7, idx8, idy8, idz8,\
er11, er12, er13, er21, er22, er23, er31, er32, er33,\
cond11, con12, cond13, cond21, con22, cond23, cond31, con32, cond33,\
Denominator,\
inverse_er11, inverse_er12, inverse_er13,\
inverse_er21, inverse_er22, inverse_er23,\
inverse_er31, inverse_er32, inverse_er33,\
cond11, cond12, cond13, cond21, cond22, cond23, cond31, cond32, cond33,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,isPML,  \
P_count) default(present) async(1)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);


		if ((i == 0 || i == nx) && periodicX_EM == false) {
			if (if_PML_Xe == false && if_PML_Xs == false) {
				continue;
			}
			else {
				idx1 = i; idx2 = i; idx3 = i; idx4 = i; idx5 = i; idx6 = i; idx7 = i; idx8 = i;	
			}
		}
		else if ((i == 0 || i == nx) && periodicX_EM == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = nx - 1; idx4 = nx - 1; idx5 = 0; idx6 = 0; idx7 = 0; idx8 = 0;
		}
		else if (i < nx + 1) {
			idx1 = i - 1; idx2 = i - 1; idx3 = i - 1; idx4 = i - 1; idx5 = i; idx6 = i; idx7 = i; idx8 = i;
		}
		if ((j == 0 || j == ny) && periodicY_EM == false) {
			if (if_PML_Ye == false && if_PML_Ys == false) {
				continue;
			}
			else {
			idy1 = j; idy2 = j; idy3 = j; idy4 = j; idy5 = j; idy6 = j; idy7 = j; idy8 = j;
			}
		}
		else if ((j == 0 || j == ny) && periodicY_EM == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0; idy5 = ny - 1; idy6 = ny - 1; idy7 = 0; idy8 = 0;
		}
		else if (j < ny + 1) {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j; idy5 = j - 1; idy6 = j - 1; idy7 = j; idy8 = j;
		}
		if ((k == 0 || k == nz) && periodicZ_EM == false) {
			if (if_PML_Ze == false && if_PML_Zs == false) {
				continue;
			}
			else {
			idz1 = k; idz2 = k; idz3 = k; idz4 = k; idz5 = k; idz6 = k; idz7 = k; idz8 = k;
			}
		}
		else if ((k == 0 || k == nz) && periodicZ_EM == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0; idz5 = nz - 1; idz6 = 0; idz7 = nz - 1; idz8 = 0;
		}
		else if (k < nz + 1) {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k; idz5 = k - 1; idz6 = k; idz7 = k - 1; idz8 = k;
		}

		dHzdy = (DHz_em_store(idx5, idy3, idz4) - DHz_em_store(idx5, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx5, idy3, idz1)) / dz;
		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.;

		dHxdz = (DHx_em_store(idx5, idy7, idz2) - DHx_em_store(idx5, idy7, idz1)) / dz;
		dHzdx = (DHz_em_store(idx5, idy7, idz2) - DHz_em_store(idx1, idy7, idz2)) / dx;
		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.;

		dHydx = (DHy_em_store(idx5, idy3, idz2) - DHy_em_store(idx1, idy3, idz2)) / dx;
		dHxdy = (DHx_em_store(idx5, idy3, idz2) - DHx_em_store(idx5, idy1, idz2)) / dy;
		Jfz = 0.; Jpz = 0.;
		dPz = 0.;

		cond11 = 0.; cond12 = 0.; cond13 = 0.; cond21 = 0.; cond22 = 0.; cond23 = 0.; cond31 = 0.; cond32 = 0.; cond33 = 0.;
		er11 = 0.; er12 = 0.; er13 = 0.; er21 = 0.; er22 = 0.; er23 = 0.; er31 = 0.; er32 = 0.; er33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;



		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			if ((idx1 < xS && true == if_PML_Xs) \
			|| (idx1 >= xE && true == if_PML_Xe) \
			|| (idy1 < yS && true == if_PML_Ys) \
			|| (idy1 >= yE && true == if_PML_Ye) \
			|| (idz1 < zS && true == if_PML_Zs) \
			|| (idz1 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx1, idy1, idz1); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx1, idy1, idz1); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx1, idy1, idz1);
		Jfy = Jfy + DJfy(idx1, idy1, idz1);
		Jfz = Jfz + DJfz(idx1, idy1, idz1);
		Jpx = Jpx + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jpy = Jpy + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jpz = Jpz + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx1, idy1, idz1);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx1, idy1, idz1);


		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			if ((idx2 < xS && true == if_PML_Xs) \
				|| (idx2 >= xE && true == if_PML_Xe) \
				|| (idy2 < yS && true == if_PML_Ys) \
				|| (idy2 >= yE && true == if_PML_Ye) \
				|| (idz2 < zS && true == if_PML_Zs) \
				|| (idz2 >= zE && true == if_PML_Ze))
			{
				isPML = true;
			}
			else
				isPML = false;

			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx2, idy2, idz2); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx2, idy2, idz2); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx2, idy2, idz2);
		Jfy = Jfy + DJfy(idx2, idy2, idz2);
		Jfz = Jfz + DJfz(idx2, idy2, idz2);
		Jpx = Jpx + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jpy = Jpy + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jpz = Jpz + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx2, idy2, idz2);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx2, idy2, idz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		////printf("mat_type3 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx3, idy3, idz3);
		if (mat_type == 0) {
			if ((idx3 < xS && true == if_PML_Xs) \
			|| (idx3 >= xE && true == if_PML_Xe) \
			|| (idy3 < yS && true == if_PML_Ys) \
			|| (idy3 >= yE && true == if_PML_Ye) \
			|| (idz3 < zS && true == if_PML_Zs) \
			|| (idz3 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx3, idy3, idz3); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx3, idy3, idz3); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx3, idy3, idz3);
		Jfy = Jfy + DJfy(idx3, idy3, idz3);
		Jfz = Jfz + DJfz(idx3, idy3, idz3);
		Jpx = Jpx + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jpy = Jpy + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jpz = Jpz + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx3, idy3, idz3);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx3, idy3, idz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		//printf("mat_type4 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx4, idy4, idz4);
		if (mat_type == 0) {
			if ((idx4 < xS && true == if_PML_Xs) \
			|| (idx4 >= xE && true == if_PML_Xe) \
			|| (idy4 < yS && true == if_PML_Ys) \
			|| (idy4 >= yE && true == if_PML_Ye) \
			|| (idz4 < zS && true == if_PML_Zs) \
			|| (idz4 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx4, idy4, idz4); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx4, idy4, idz4); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx4, idy4, idz4);
		Jfy = Jfy + DJfy(idx4, idy4, idz4);
		Jfz = Jfz + DJfz(idx4, idy4, idz4);
		Jpx = Jpx + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jpy = Jpy + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jpz = Jpz + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx4, idy4, idz4);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx4, idy4, idz4);


		//----------Polarization 5------//
		mat_type = pt_glb->material_cell(idx5, idy5, idz5);
		if (mat_type == 0) {
			if ((idx5 < xS && true == if_PML_Xs) \
			|| (idx5 >= xE && true == if_PML_Xe) \
			|| (idy5 < yS && true == if_PML_Ys) \
			|| (idy5 >= yE && true == if_PML_Ye) \
			|| (idz5 < zS && true == if_PML_Zs) \
			|| (idz5 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx5, idy5, idz5); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx5, idy5, idz5); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx5, idy5, idz5); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx5, idy5, idz5);
		Jfy = Jfy + DJfy(idx5, idy5, idz5);
		Jfz = Jfz + DJfz(idx5, idy5, idz5);
		Jpx = Jpx + Jpx_n1_store(idx5, idy5, idz5) + Jpx_n2_store(idx5, idy5, idz5) + Jpx_n3_store(idx5, idy5, idz5);
		Jpy = Jpy + Jpy_n1_store(idx5, idy5, idz5) + Jpy_n2_store(idx5, idy5, idz5) + Jpy_n3_store(idx5, idy5, idz5);
		Jpz = Jpz + Jpz_n1_store(idx5, idy5, idz5) + Jpz_n2_store(idx5, idy5, idz5) + Jpz_n3_store(idx5, idy5, idz5);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx5, idy5, idz5);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx5, idy5, idz5);



		//----------Polarization 6------//
		mat_type = pt_glb->material_cell(idx6, idy6, idz6);
		//printf("mat_type6 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx6, idy6, idz6);

		if (mat_type == 0) {
			if ((idx6 < xS && true == if_PML_Xs) \
			|| (idx6 >= xE && true == if_PML_Xe) \
			|| (idy6 < yS && true == if_PML_Ys) \
			|| (idy6 >= yE && true == if_PML_Ye) \
			|| (idz6 < zS && true == if_PML_Zs) \
			|| (idz6 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx6, idy6, idz6); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx6, idy6, idz6); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx6, idy6, idz6); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx6, idy6, idz6);
		Jfy = Jfy + DJfy(idx6, idy6, idz6);
		Jfz = Jfz + DJfz(idx6, idy6, idz6);
		Jpx = Jpx + Jpx_n1_store(idx6, idy6, idz6) + Jpx_n2_store(idx6, idy6, idz6) + Jpx_n3_store(idx6, idy6, idz6);
		Jpy = Jpy + Jpy_n1_store(idx6, idy6, idz6) + Jpy_n2_store(idx6, idy6, idz6) + Jpy_n3_store(idx6, idy6, idz6);
		Jpz = Jpz + Jpz_n1_store(idx6, idy6, idz6) + Jpz_n2_store(idx6, idy6, idz6) + Jpz_n3_store(idx6, idy6, idz6);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx6, idy6, idz6);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx6, idy6, idz6);



		//----------Polarization 7------//
		mat_type = pt_glb->material_cell(idx7, idy7, idz7);
		//printf("mat_type7 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx7, idy7, idz7);

		if (mat_type == 0) {
			if ((idx7 < xS && true == if_PML_Xs) \
			|| (idx7 >= xE && true == if_PML_Xe) \
			|| (idy7 < yS && true == if_PML_Ys) \
			|| (idy7 >= yE && true == if_PML_Ye) \
			|| (idz7 < zS && true == if_PML_Zs) \
			|| (idz7 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx7, idy7, idz7); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx7, idy7, idz7); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx7, idy7, idz7); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx7, idy7, idz7);
		Jfy = Jfy + DJfy(idx7, idy7, idz7);
		Jfz = Jfz + DJfz(idx7, idy7, idz7);
		Jpx = Jpx + Jpx_n1_store(idx7, idy7, idz7) + Jpx_n2_store(idx7, idy7, idz7) + Jpx_n3_store(idx7, idy7, idz7);
		Jpy = Jpy + Jpy_n1_store(idx7, idy7, idz7) + Jpy_n2_store(idx7, idy7, idz7) + Jpy_n3_store(idx7, idy7, idz7);
		Jpz = Jpz + Jpz_n1_store(idx7, idy7, idz7) + Jpz_n2_store(idx7, idy7, idz7) + Jpz_n3_store(idx7, idy7, idz7);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx7, idy7, idz7);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx7, idy7, idz7);


		//----------Polarization 8------//
		mat_type = pt_glb->material_cell(idx8, idy8, idz8);
		//printf("mat_type8 = %d\t$%ld\t%ld\t%ld\n", mat_type, idx8, idy8, idz8);
		
		if (mat_type == 0) {
			if ((idx8 < xS && true == if_PML_Xs) \
			|| (idx8 >= xE && true == if_PML_Xe) \
			|| (idy8 < yS && true == if_PML_Ys) \
			|| (idy8 >= yE && true == if_PML_Ye) \
			|| (idz8 < zS && true == if_PML_Zs) \
			|| (idz8 >= zE && true == if_PML_Ze)) 
			{
				isPML = true;
			}
			else
				isPML = false;
			
			if (isPML) {
				er11 += pt_glb->PML_er11; er22 += pt_glb->PML_er22; er33 += pt_glb->PML_er33;
			}
			else {
				er11 = er11 + 1.; 	er22 = er22 + 1.; 	er33 = er33 + 1.;
			}
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(idx8, idy8, idz8); //RK
				dPy = dPy + pt_fe->dpy_glb_rk4(idx8, idy8, idz8); //RK
				dPz = dPz + pt_fe->dpz_glb_rk4(idx8, idy8, idz8); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
			cond11 = cond11 + mat->conductivity11; cond12 = cond12 + mat->conductivity12; cond13 = cond13 + mat->conductivity13;
			cond21 = cond21 + mat->conductivity21; cond22 = cond22 + mat->conductivity22; cond23 = cond23 + mat->conductivity23;
			cond31 = cond31 + mat->conductivity31; cond32 = cond32 + mat->conductivity32; cond33 = cond33 + mat->conductivity33;
		}
		Jfx = Jfx + DJfx(idx8, idy8, idz8);
		Jfy = Jfy + DJfy(idx8, idy8, idz8);
		Jfz = Jfz + DJfz(idx8, idy8, idz8);
		Jpx = Jpx + Jpx_n1_store(idx8, idy8, idz8) + Jpx_n2_store(idx8, idy8, idz8) + Jpx_n3_store(idx8, idy8, idz8);
		Jpy = Jpy + Jpy_n1_store(idx8, idy8, idz8) + Jpy_n2_store(idx8, idy8, idz8) + Jpy_n3_store(idx8, idy8, idz8);
		Jpz = Jpz + Jpz_n1_store(idx8, idy8, idz8) + Jpz_n2_store(idx8, idy8, idz8) + Jpz_n3_store(idx8, idy8, idz8);
		Jishex = Jishex + pt_mag->Jx_ISHE(idx8, idy8, idz8);
		Jishey = Jishey + pt_mag->Jx_ISHE(idx8, idy8, idz8);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jfx = Jfx / 8.;	Jfy = Jfy / 8.;	Jfz = Jfz / 8.;
		Jpx = Jpx / 8.;	Jpy = Jpy / 8.;	Jpz = Jpz / 8.;
		Jishex = Jishex / 8.;	Jishey = Jishey / 8.;
		dPx = dPx / P_count;	dPy = dPy / P_count;	dPz = dPz / P_count;

		cond11 = cond11 / 8.; cond12 = cond12 / 8.; cond13 = cond13 / 8.;
		cond21 = cond21 / 8.; cond22 = cond22 / 8.; cond23 = cond23 / 8.;
		cond31 = cond31 / 8.; cond32 = cond32 / 8.; cond33 = cond33 / 8.;

		er11 = er11 / 8.; er12 = er12 / 8.; er13 = er13 / 8.;
		er21 = er21 / 8.; er22 = er22 / 8.; er23 = er23 / 8.;
		er31 = er31 / 8.; er32 = er32 / 8.; er33 = er33 / 8.;


		Denominator = -er13 * er22 * er31 + er12 * er23 * er31 + er13 * er21 * er32 - \
			er11 * er23 * er32 - er12 * er21 * er33 + er11 * er22 * er33;
		inverse_er11 = (-er23 * er32 + er22 * er33) / Denominator;
		inverse_er12 = (er13 * er32 - er12 * er33) / Denominator;
		inverse_er13 = (-er13 * er22 + er12 * er23) / Denominator;
		inverse_er21 = (er23 * er31 - er21 * er33) / Denominator;
		inverse_er22 = (-er13 * er31 + er11 * er33) / Denominator;
		inverse_er23 = (er13 * er21 - er11 * er23) / Denominator;
		inverse_er31 = (-er22 * er31 + er21 * er32) / Denominator;
		inverse_er32 = (er12 * er31 - er11 * er32) / Denominator;
		inverse_er33 = (-er12 * er21 + er11 * er22) / Denominator;

		isPML = false;

		if (((true == if_PML_Xs && i < xS) || (true == if_PML_Xe && i >= xE)) \
			|| ((true == if_PML_Ys && j < yS) || (true == if_PML_Ye && j >= yE)) \
			|| ((true == if_PML_Zs && k < zS) || (true == if_PML_Ze && k >= zE)))
		{
			isPML = true;

			if (i < nx)
				DEx_em_t4(idx8, j, k) = pt_glb->dt / kappa_y_np1(0, j, 0) \
					* (dHzdy - dHydz - (sigma_y_np1(0, j, 0) / e0) * Dx_PML_store(idx8, j, k) \
					- Jfx - Jpx - Jishex - dPx / pt_glb->dt);

			if (j < ny)
				DEy_em_t4(i, idy8, k) = pt_glb->dt / kappa_z_np1(0, 0, k) \
					* (dHxdz - dHzdx - (sigma_z_np1(0, 0, k) / e0) * Dy_PML_store(i, idy8, k) \
						- Jfy - Jpy - Jishey - dPy / pt_glb->dt);

			if (k < nz)
				DEz_em_t4(i, j, idz8) = pt_glb->dt / kappa_x_np1(i, 0, 0) \
					* (dHydx - dHxdy - (sigma_x_np1(i, 0, 0) / e0) * Dz_PML_store(i, j, idz8) \
					- Jfz - Jpz - dPz / pt_glb->dt);
		}

		if (true == isPML)
		{
			if (i < nx) {
				dDEx_em_rk4(idx8, j, k) = (pt_glb->dt / kappa_z_np1(0, 0, k) / e0) \
					* (inverse_er11*(kappa_x_n(idx8, 0, 0) * DEx_em_t4(idx8, j, k) / pt_glb->dt \
						+ sigma_x_n(idx8, 0, 0) * Dx_PML_store(idx8, j, k) / e0)\
						- (sigma_z_np1(0, 0, k) * DEx_em_store(idx8, j, k)));
			}

			if (j < ny) {
				dDEy_em_rk4(i, idy8, k) = (pt_glb->dt / kappa_x_np1(i, 0, 0) / e0) \
					* (inverse_er22*(kappa_y_n(0, idy8, 0) * DEy_em_t4(i, idy8, k) / pt_glb->dt \
						+ sigma_y_n(0, idy8, 0) * Dy_PML_store(i, idy8, k) / e0)\
						- (sigma_x_np1(i, 0, 0) * DEy_em_store(i, idy8, k)));
			}

			if (k < nz) {
				dDEz_em_rk4(i, j, idz8) = (pt_glb->dt / kappa_y_np1(0, j, 0) / e0) \
					* (inverse_er33*(kappa_z_n(0, 0, idz8) * DEz_em_t4(i, j, idz8) / pt_glb->dt \
						+ sigma_z_n(0, 0, idz8) * Dz_PML_store(i, j, idz8) / e0) \
						- (sigma_y_np1(0, j, 0) * DEz_em_store(i, j, idz8)));
			}
		}
		else
		{
			if (i < nx) {
				dDEx_em_rk4(idx8, j, k) = \
					pt_glb->dt / e0 * inverse_er11 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er12 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er13 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (j < ny) {
				dDEy_em_rk4(i, idy8, k) = \
					pt_glb->dt / e0 * inverse_er21 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er22 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er23 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}

			if (k < nz) {
				dDEz_em_rk4(i, j, idz8) = \
					pt_glb->dt / e0 * inverse_er31 * \
					(dHzdy - dHydz - Jfx - Jpx - Jishex - \
						cond11 * DEx_em_store(idx8, j, k) - cond12 * DEy_em_store(i, idy8, k) - cond13 * DEz_em_store(i, j, idz8) - \
						dPx / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er32 * \
					(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
						cond21 * DEx_em_store(idx8, j, k) - cond22 * DEy_em_store(i, idy8, k) - cond23 * DEz_em_store(i, j, idz8) - \
						dPy / pt_glb->dt) + \
					pt_glb->dt / e0 * inverse_er33 * \
					(dHydx - dHxdy - Jfz - Jpz - \
						cond31 * DEx_em_store(idx8, j, k) - cond32 * DEy_em_store(i, idy8, k) - cond33 * DEz_em_store(i, j, idz8) - \
						dPz / pt_glb->dt); //RK
			}
		}
	}
	//-----------END X Y Z component-----------//
}

void EMdynamic_system::get_dH_RK1() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	bool isPML = false;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count, isPML) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);
		
		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size + 1 && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		if (i == 0) {
			if (periodicX_EM == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (periodicX_EM == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk1(idx1, idy1, idz1) * Ms1 \
					+ pt_mag->dmx_AFM1_glb_rk1(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk1(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk1(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk1(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk1(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (j < ny && k < nz) {
				BHx_em_t1(i, j, k) = -pt_glb->dt / kappa_y_n(0, j, 0) \
					* (dEzdy - dEydz + (sigma_y_n(0, j, 0) / e0) * Bx_PML_store(i, j, k)) \
						- (dMx1 + dMx2) / M_count;

				dDHx_em_rk1(i, j, k) = (pt_glb->dt / kappa_z_n(0, 0, k)) \
					* ((kappa_x_np1(i, 0, 0) * BHx_em_t1(i, j, k) / pt_glb->dt \
						+ sigma_x_np1(i, 0, 0) * Bx_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_z_n(0, 0, k) * DHx_em_store(i, j, k) / e0));
			}
		}
		else
			dDHx_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		
			//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count,isPML) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size + 1 && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		idx1 = i; idx2 = i;
		if (j == 0) {
			if (periodicY_EM == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (periodicY_EM == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk1(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk1(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk1(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk1(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk1(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk1(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (i < nx && k < nz) {
				BHy_em_t1(i, j, k) = -pt_glb->dt / kappa_z_n(0, 0, k) \
					* (dExdz - dEzdx + (sigma_z_n(0, 0, k) / e0) * By_PML_store(i, j, k)) \
						- (dMy1 + dMy2) / M_count;

				dDHy_em_rk1(i, j, k) = (pt_glb->dt / kappa_x_n(i, 0, 0)) \
					* ((kappa_y_np1(0, j, 0) * BHy_em_t1(i, j, k) / pt_glb->dt \
						+ sigma_y_np1(0, j, 0) * By_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_x_n(i, 0, 0) * DHy_em_store(i, j, k) / e0));
			}
		}
		else
			dDHy_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count,isPML) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= nz - PML_size + 1 && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (periodicZ_EM == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk1(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk1(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk1(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (periodicZ_EM == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk1(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk1(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk1(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (i < nx && j < ny) {
				BHz_em_t1(i, j, k) = -pt_glb->dt / kappa_x_n(i, 0, 0) \
					* (dEydx - dExdy + (sigma_x_n(i, 0, 0) / e0) * Bz_PML_store(i, j, k)) \
						- (dMz1 + dMz2) / M_count;

				dDHz_em_rk1(i, j, k) = (pt_glb->dt / kappa_y_n(0, j, 0)) \
					* ((kappa_z_np1(0, 0, k) * BHz_em_t1(i, j, k) / pt_glb->dt \
						+ sigma_z_np1(0, 0, k) * Bz_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_y_n(0, j, 0) * DHz_em_store(i, j, k) / e0));
			}
		}
		else
			dDHz_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}

void EMdynamic_system::get_dH_RK2() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	bool isPML;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count,isPML) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size + 1 && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		if (i == 0) {
			if (periodicX_EM == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (periodicX_EM == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk2(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmx_AFM1_glb_rk2(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk2(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk2(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk2(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk2(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (j < ny && k < nz) {
				BHx_em_t2(i, j, k) = -pt_glb->dt / kappa_y_n(0, j, 0) \
					* (dEzdy - dEydz + (sigma_y_n(0, j, 0) / e0) * Bx_PML_store(i, j, k)) \
						- (dMx1 + dMx2) / M_count;

				dDHx_em_rk2(i, j, k) = (pt_glb->dt / kappa_z_n(0, 0, k)) \
					* ((kappa_x_np1(i, 0, 0) * BHx_em_t2(i, j, k) / pt_glb->dt \
						+ sigma_x_np1(i, 0, 0) * Bx_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_z_n(0, 0, k) * DHx_em_store(i, j, k) / e0));
			}
		}
		else
			dDHx_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count,isPML) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size + 1 && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		idx1 = i; idx2 = i;
		if (j == 0) {
			if (periodicY_EM == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (periodicY_EM == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk2(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk2(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk2(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk2(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk2(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk2(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		
		if (true == if_PML && true == isPML) {
			if (i < nx && k < nz) {
				BHy_em_t2(i, j, k) = -pt_glb->dt / kappa_z_n(0, 0, k) \
					* (dExdz - dEzdx + (sigma_z_n(0, 0, k) / e0) * By_PML_store(i, j, k)) \
						- (dMy1 + dMy2) / M_count;

				dDHy_em_rk2(i, j, k) = (pt_glb->dt / kappa_x_n(i, 0, 0)) \
					* ((kappa_y_np1(0, j, 0) * BHy_em_t2(i, j, k) / pt_glb->dt \
						+ sigma_y_np1(0, j, 0) * By_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_x_n(i, 0, 0) * DHy_em_store(i, j, k) / e0));
			}
		}
		else
			dDHy_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count,isPML) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= nz - PML_size + 1 && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (periodicZ_EM == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk2(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk2(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk2(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (periodicZ_EM == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk2(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk2(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk2(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (i < nx && j < ny) {
				BHz_em_t2(i, j, k) = -pt_glb->dt / kappa_x_n(i, 0, 0) \
					* (dEydx - dExdy + (sigma_x_n(i, 0, 0) / e0) * Bz_PML_store(i, j, k)) \
						- (dMz1 + dMz2) / M_count;

				dDHz_em_rk2(i, j, k) = (pt_glb->dt / kappa_y_n(0, j, 0)) \
					* ((kappa_z_np1(0, 0, k) * BHz_em_t2(i, j, k) / pt_glb->dt \
						+ sigma_z_np1(0, 0, k) * Bz_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_y_n(0, j, 0) * DHz_em_store(i, j, k) / e0));
			}
		}
		else
			dDHz_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}

void EMdynamic_system::get_dH_RK3() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	bool isPML;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count,isPML) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size + 1 && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		if (i == 0) {
			if (periodicX_EM == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (periodicX_EM == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk3(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmx_AFM1_glb_rk3(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk3(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk3(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk3(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk3(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (j < ny && k < nz) {
				BHx_em_t3(i, j, k) = -pt_glb->dt / kappa_y_n(0, j, 0) \
					* (dEzdy - dEydz + (sigma_y_n(0, j, 0) / e0) * Bx_PML_store(i, j, k)) \
						- (dMx1 + dMx2) / M_count;

				dDHx_em_rk3(i, j, k) = (pt_glb->dt / kappa_z_n(0, 0, k)) \
					* ((kappa_x_np1(i, 0, 0) * BHx_em_t3(i, j, k) / pt_glb->dt \
						+ sigma_x_np1(i, 0, 0) * Bx_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_z_n(0, 0, k) * DHx_em_store(i, j, k) / e0));
			}
		}
		else
			dDHx_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK

		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count,isPML) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size + 1 && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		idx1 = i; idx2 = i;
		if (j == 0) {
			if (periodicY_EM == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (periodicY_EM == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk3(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk3(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk3(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk3(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk3(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk3(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (i < nx && k < nz) {
				BHy_em_t3(i, j, k) = -pt_glb->dt / kappa_z_n(0, 0, k) \
					* (dExdz - dEzdx + (sigma_z_n(0, 0, k) / e0) * By_PML_store(i, j, k)) \
						- (dMy1 + dMy2) / M_count;

				dDHy_em_rk3(i, j, k) = (pt_glb->dt / kappa_x_n(i, 0, 0)) \
					* ((kappa_y_np1(0, j, 0) * BHy_em_t3(i, j, k) / pt_glb->dt \
						+ sigma_y_np1(0, j, 0) * By_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_x_n(i, 0, 0) * DHy_em_store(i, j, k) / e0));
			}
		}
		else
			dDHy_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count,isPML) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= nz - PML_size + 1 && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (periodicZ_EM == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk3(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk3(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk3(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (periodicZ_EM == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk3(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk3(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk3(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (i < nx && j < ny) {
				BHz_em_t3(i, j, k) = -pt_glb->dt / kappa_x_n(i, 0, 0) \
					* (dEydx - dExdy + (sigma_x_n(i, 0, 0) / e0) * Bz_PML_store(i, j, k)) \
						- (dMz1 + dMz2) / M_count;

				dDHz_em_rk3(i, j, k) = (pt_glb->dt / kappa_y_n(0, j, 0)) \
					* ((kappa_z_np1(0, 0, k) * BHz_em_t3(i, j, k) / pt_glb->dt \
						+ sigma_z_np1(0, 0, k) * Bz_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_y_n(0, j, 0) * DHz_em_store(i, j, k) / e0));
			}
		}
		else
			dDHz_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}

void EMdynamic_system::get_dH_RK4() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	bool isPML;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count, isPML) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size + 1 && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		if (i == 0) {
			if (periodicX_EM == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (periodicX_EM == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk4(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmx_AFM1_glb_rk4(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk4(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk4(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk4(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk4(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (j < ny && k < nz) {
				BHx_em_t4(i, j, k) = -pt_glb->dt / kappa_y_n(0, j, 0) \
					* (dEzdy - dEydz + (sigma_y_n(0, j, 0) / e0) * Bx_PML_store(i, j, k)) \
						- (dMx1 + dMx2) / M_count;

				dDHx_em_rk4(i, j, k) = (pt_glb->dt / kappa_z_n(0, 0, k)) \
					* ((kappa_x_np1(i, 0, 0) * BHx_em_t4(i, j, k) / pt_glb->dt \
						+ sigma_x_np1(i, 0, 0) * Bx_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_z_n(0, 0, k) * DHx_em_store(i, j, k) / e0));
			}
		}
		else
			dDHx_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count,isPML) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size + 1 && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= zE && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (periodicY_EM == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (periodicY_EM == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk4(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk4(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk4(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk4(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk4(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk4(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (i < nx && k < nz) {
				BHy_em_t4(i, j, k) = -pt_glb->dt / kappa_z_n(0, 0, k) \
					* (dExdz - dEzdx + (sigma_z_n(0, 0, k) / e0) * By_PML_store(i, j, k)) \
						- (dMy1 + dMy2) / M_count;

				dDHy_em_rk4(i, j, k) = (pt_glb->dt / kappa_x_n(i, 0, 0)) \
					* ((kappa_y_np1(0, j, 0) * BHy_em_t4(i, j, k) / pt_glb->dt \
						+ sigma_y_np1(0, j, 0) * By_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_x_n(i, 0, 0) * DHy_em_store(i, j, k) / e0));
			}
		}
		else
			dDHy_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count,isPML) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i < PML_size && true == if_PML_Xs) \
			|| (i >= nx - PML_size && true == if_PML_Xe) \
			|| (j < PML_size && true == if_PML_Ys) \
			|| (j >= ny - PML_size && true == if_PML_Ye) \
			|| (k < PML_size && true == if_PML_Zs) \
			|| (k >= nz - PML_size + 1 && true == if_PML_Ze)) 
		{
			isPML = true;
		}
		else
			isPML = false;
		
		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (periodicZ_EM == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk4(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk4(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk4(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (periodicZ_EM == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk4(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk4(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk4(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}

		if (true == if_PML && true == isPML) {
			if (i < nx && j < ny) {
				BHz_em_t4(i, j, k) = -pt_glb->dt / kappa_x_n(i, 0, 0) \
					* (dEydx - dExdy + (sigma_x_n(i, 0, 0) / e0) * Bz_PML_store(i, j, k)) \
						- (dMz1 + dMz2) / M_count;

				dDHz_em_rk4(i, j, k) = (pt_glb->dt / kappa_y_n(0, j, 0)) \
					* ((kappa_z_np1(0, 0, k) * BHz_em_t4(i, j, k) / pt_glb->dt \
						+ sigma_z_np1(0, 0, k) * Bz_PML_store(i, j, k) / e0) / mu0 \
						- (sigma_y_n(0, j, 0) * DHz_em_store(i, j, k) / e0));
			}
		}
		else
			dDHz_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}


void EMdynamic_system::update_DH_RK1()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em_store(id) = DHx_em(id) + dDHx_em_rk1(id) * 0.5;
			if (true == if_PML)
				Bx_PML_store(id) = Bx_PML(id) + BHx_em_t1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk1(id) * 0.5;
			if (true == if_PML)
				By_PML_store(id) = By_PML(id) + BHy_em_t1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk1(id) * 0.5;
			if (true == if_PML)
				Bz_PML_store(id) = Bz_PML(id) + BHz_em_t1(id) * 0.5;
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DH_RK2()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em_store(id) = DHx_em(id) + dDHx_em_rk2(id) * 0.5;
			if (true == if_PML)
				Bx_PML_store(id) = Bx_PML(id) + BHx_em_t2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk2(id) * 0.5;
			if (true == if_PML)
				By_PML_store(id) = By_PML(id) + BHy_em_t2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk2(id) * 0.5;
			if (true == if_PML)
				Bz_PML_store(id) = Bz_PML(id) + BHz_em_t2(id) * 0.5;
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DH_RK3()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em_store(id) = DHx_em(id) + dDHx_em_rk3(id);
			if (true == if_PML)
				Bx_PML_store(id) = Bx_PML(id) + BHx_em_t3(id) ;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk3(id);
			if (true == if_PML)
				By_PML_store(id) = By_PML(id) + BHy_em_t3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk3(id);
			if (true == if_PML)
				Bz_PML_store(id) = Bz_PML(id) + BHz_em_t3(id);
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DH()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em(id) = DHx_em(id) + dDHx_em_rk1(id) / 6. + dDHx_em_rk2(id) / 3. + dDHx_em_rk3(id) / 3. + dDHx_em_rk4(id) / 6.;
			DHx_em_store(id) = DHx_em(id);
			if (true == if_PML)
			{
				Bx_PML(id) = Bx_PML(id) + BHx_em_t1(id) / 6. + BHx_em_t2(id) / 3. + BHx_em_t3(id) / 3. + BHx_em_t4(id) / 6.;
				Bx_PML_store(id) = Bx_PML(id);
			}
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em(id) = DHy_em(id) + dDHy_em_rk1(id) / 6. + dDHy_em_rk2(id) / 3. + dDHy_em_rk3(id) / 3. + dDHy_em_rk4(id) / 6.;
			DHy_em_store(id) = DHy_em(id);
			if (true == if_PML)
			{
				By_PML(id) = By_PML(id) + BHy_em_t1(id) / 6. + BHy_em_t2(id) / 3. + BHy_em_t3(id) / 3. + BHy_em_t4(id) / 6.;
				By_PML_store(id) = By_PML(id);
			}
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em(id) = DHz_em(id) + dDHz_em_rk1(id) / 6. + dDHz_em_rk2(id) / 3. + dDHz_em_rk3(id) / 3. + dDHz_em_rk4(id) / 6.;
			DHz_em_store(id) = DHz_em(id);
			if (true == if_PML)
			{
				Bz_PML(id) = Bz_PML(id) + BHz_em_t1(id) / 6. + BHz_em_t2(id) / 3. + BHz_em_t3(id) / 3. + BHz_em_t4(id) / 6.;
				Bz_PML_store(id) = Bz_PML(id);
			}
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DE_RK1()
{
#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em_store(id) = DEx_em(id) + dDEx_em_rk1(id) * 0.5;
			//printf("DEx1: %le\t%le\t%ld\n", DEx_em(id), dDEx_em_rk1(id), id);
			if (true == if_PML)
				Dx_PML_store(id) = Dx_PML(id) + DEx_em_t1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk1(id) * 0.5;
			//printf("DEy1: %le\t%le\t%ld\n", DEy_em(id), dDEy_em_rk1(id), id);
			if (true == if_PML)
				Dy_PML_store(id) = Dy_PML(id) + DEy_em_t1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk1(id) * 0.5;
			//printf("DEz1: %le\t%le\t%ld\n", DEz_em(id), dDEz_em_rk1(id), id);
			if (true == if_PML)
				Dz_PML_store(id) = Dz_PML(id) + DEz_em_t1(id) * 0.5;
		}
	}

if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
	{
		update_planeEM_half();
	}
}

if (pt_glb->if_periodic_allsurface == false && if_PML == false) {
	update_DE_Boundary_half();
}

#pragma acc parallel default(present) async(8)
{
	update_DE_cell();
}


#pragma acc parallel default(present) async(9)
{
	update_Jp_RK1();
}
}

void EMdynamic_system::update_DE_RK2()
{
#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em_store(id) = DEx_em(id) + dDEx_em_rk2(id) * 0.5;
			if (true == if_PML)
				Dx_PML_store(id) = Dx_PML(id) + DEx_em_t2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk2(id) * 0.5;
			if (true == if_PML)
				Dy_PML_store(id) = Dy_PML(id) + DEy_em_t2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk2(id) * 0.5;
			if (true == if_PML)
				Dz_PML_store(id) = Dz_PML(id) + DEz_em_t2(id) * 0.5;
		}
	}

if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
	{
		update_planeEM_half();
	}
}

if (pt_glb->if_periodic_allsurface == false && if_PML == false) {
	update_DE_Boundary_half();
}

#pragma acc parallel default(present) async(8)
{
	update_DE_cell();
}

#pragma acc parallel default(present) async(9)
{
	update_Jp_RK2();
}
}

void EMdynamic_system::update_DE_RK3()
{
#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em_store(id) = DEx_em(id) + dDEx_em_rk3(id);
			if (true == if_PML)
				Dx_PML_store(id) = Dx_PML(id) + DEx_em_t3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk3(id);
			if (true == if_PML)
				Dy_PML_store(id) = Dy_PML(id) + DEy_em_t3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk3(id);
			if (true == if_PML)
				Dz_PML_store(id) = Dz_PML(id) + DEz_em_t3(id);
		}
	}

if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
	{
		update_planeEM_full();
	}
}

if (pt_glb->if_periodic_allsurface == false && if_PML == false) {
	update_DE_Boundary_full();
}

#pragma acc parallel default(present) async(8)
{
	update_DE_cell();
}

#pragma acc parallel default(present) async(9)
{
	update_Jp_RK3();
}
}

void EMdynamic_system::update_DE()
{
	if (pt_glb->if_periodic_allsurface == false && if_PML == false) {
		transfer_pointer();
	}

#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em(id) = DEx_em(id) + dDEx_em_rk1(id) / 6. + dDEx_em_rk2(id) / 3. + dDEx_em_rk3(id) / 3. + dDEx_em_rk4(id) / 6.;
			if (true == if_PML)
				Dx_PML(id) = Dx_PML(id) + DEx_em_t1(id) / 6. + DEx_em_t2(id) / 3. + DEx_em_t3(id) / 3. + DEx_em_t4(id) / 6.;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em(id) = DEy_em(id) + dDEy_em_rk1(id) / 6. + dDEy_em_rk2(id) / 3. + dDEy_em_rk3(id) / 3. + dDEy_em_rk4(id) / 6.;
			if (true == if_PML)
				Dy_PML(id) = Dy_PML(id) + DEy_em_t1(id) / 6. + DEy_em_t2(id) / 3. + DEy_em_t3(id) / 3. + DEy_em_t4(id) / 6.;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em(id) = DEz_em(id) + dDEz_em_rk1(id) / 6. + dDEz_em_rk2(id) / 3. + dDEz_em_rk3(id) / 3. + dDEz_em_rk4(id) / 6.;
			if (true == if_PML)
				Dz_PML(id) = Dz_PML(id) + DEz_em_t1(id) / 6. + DEz_em_t2(id) / 3. + DEz_em_t3(id) / 3. + DEz_em_t4(id) / 6.;
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM();
		}
	}

	if (pt_glb->if_periodic_allsurface == false && if_PML == false) {
		update_DE_Boundary();
	}

#pragma acc parallel default(present) async(8)
	{
		DEx_em_store = DEx_em;
		DEy_em_store = DEy_em;
		DEz_em_store = DEz_em;

		if (true == if_PML) {
			Dx_PML_store = Dx_PML;
			Dy_PML_store = Dy_PML;
			Dz_PML_store = Dz_PML;
		}
	}

#pragma acc parallel default(present) async(8)
	{
		update_DE_cell();
	}

#pragma acc parallel default(present) async(9)
	{
		update_Jp();
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_DH_cell() {
	long int i, j, k;
	unsigned int mat_type;
	material* mat;

#pragma acc loop gang vector private(i,j,k,mat_type,mat)
	for (long int id = 0; id < n; id++) {
		i = id / (ny * nz);
		j = (id - i * (ny * nz)) / nz;
		k = id - i * (ny * nz) - j * nz;

		DHx_em_cell(i, j, k) = (DHx_em_store(i, j, k) + DHx_em_store(i + 1, j, k)) / 2.;
		DHy_em_cell(i, j, k) = (DHy_em_store(i, j, k) + DHy_em_store(i, j + 1, k)) / 2.;
		if (pt_glb->if_mag_1Dmodel == true && pt_glb->if_EM_fromM == true) {
			mat_type = pt_glb->material_cell(id);
			if (mat_type == 0) {
				DHz_em_cell(id) = 0.;
			}
			else {
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				DHz_em_cell(id) = -1. * mat->Ms * pt_mag->mz_glb_store(id) \
					- 1. * mat->Ms_AFM1 * pt_mag->mz_AFM1_glb_store(id) \
					- 1. * mat->Ms_AFM2 * pt_mag->mz_AFM2_glb_store(id) \
					- pt_mag->Hz_stat(id);
			}
		}
		else {
			DHz_em_cell(i, j, k) = (DHz_em_store(i, j, k) + DHz_em_store(i, j, k + 1)) / 2.;
		}
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_DE_cell() {
	long int i, j, k;

#pragma acc loop gang vector private(i,j,k)
	for (long int id = 0; id < n; id++) {
		i = id / (ny * nz);
		j = (id - i * (ny * nz)) / nz;
		k = id - i * (ny * nz) - j * nz;

		DEx_em_cell(i, j, k) = (DEx_em_store(i, j, k) + DEx_em_store(i, j + 1, k) + DEx_em_store(i, j, k + 1) + DEx_em_store(i, j + 1, k + 1)) / 4.;
		DEy_em_cell(i, j, k) = (DEy_em_store(i, j, k) + DEy_em_store(i + 1, j, k) + DEy_em_store(i, j, k + 1) + DEy_em_store(i + 1, j, k + 1)) / 4.;
		DEz_em_cell(i, j, k) = (DEz_em_store(i, j, k) + DEz_em_store(i + 1, j, k) + DEz_em_store(i, j + 1, k) + DEz_em_store(i + 1, j + 1, k)) / 4.;
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_planeEM_half() {
	unsigned long i, j;

#pragma acc loop gang vector private(i,j)
	for (unsigned long id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		DEx_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[0] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[0] * (pt_glb->time_device - 0.5 * pt_glb->dt));

		DEy_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[1] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[1] * (pt_glb->time_device - 0.5 * pt_glb->dt));

		DEz_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[2] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[2] * (pt_glb->time_device - 0.5 * pt_glb->dt));
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_planeEM_full() {
	unsigned long i, j;

#pragma acc loop gang vector private(i,j)
	for (unsigned long id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		DEx_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[0] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[0] * pt_glb->time_device);

		DEy_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[1] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[1] * pt_glb->time_device);

		DEz_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[2] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[2] * pt_glb->time_device);
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_planeEM() {
	unsigned long i, j;

#pragma acc loop gang vector private(i,j)
	for (unsigned long id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		DEx_em(i, j, planeEM_source_z) = pt_glb->planeEM_E[0] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[0] * pt_glb->time_device);

		DEy_em(i, j, planeEM_source_z) = pt_glb->planeEM_E[1] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[1] * pt_glb->time_device);

		DEz_em(i, j, planeEM_source_z) = pt_glb->planeEM_E[2] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[2] * pt_glb->time_device);
	}
}

void EMdynamic_system::update_Jf_input() {
	unsigned long int i, j, k;
	double temporal_var, spatial_var;

#pragma acc parallel default(present) async(7)
	{
#pragma acc loop gang vector private(i,j,k, temporal_var, spatial_var)
		for (unsigned long id = 0; id < Jf_n; id++) {
			i = id / (Jf_ny * Jf_nz);
			j = (id - i * (Jf_ny * Jf_nz)) / Jf_nz;
			k = id - i * (Jf_ny * Jf_nz) - j * Jf_nz;

			i = i + pt_glb->Jfin_xi;
			j = j + pt_glb->Jfin_yi;
			k = k + pt_glb->Jfin_zi;

			if (pt_glb->Jf_input_type < 3) {
				temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - pt_glb->dt));
				spatial_var = PI / (static_cast<double>(Jf_nz) * dz) * (static_cast<double>(k - pt_glb->Jfin_zi) + 0.5) * dz;
				if (pt_glb->Jf_input_type == 0) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
				}
				else if (pt_glb->Jf_input_type == 1) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
				}
				else if (pt_glb->Jf_input_type == 2) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
				}
			}
			else if (pt_glb->Jf_input_type == 3) {
				if ((pt_glb->time_device - pt_glb->dt) <= static_cast<double>(pt_glb->num_Jf_cycles) / pt_glb->Jf_input_freq) {
					temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - pt_glb->dt));
				}
				else {
					temporal_var = 0.;
				}

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
			}
			else if (pt_glb->Jf_input_type == 4) {
				temporal_var = pt_glb->time_device - pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 5) {
				temporal_var = pt_glb->time_device - pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 6) {
				temporal_var = pt_glb->time_device - pt_glb->dt;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = cos(PI * pt_glb->Jf_rotate_xy / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * pt_glb->Jf_rotate_xy / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfx(i, j, k) = cos(PI * (pt_glb->Jf_rotate_xy + 90.) / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * (pt_glb->Jf_rotate_xy + 90.) / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
			}
			else if (pt_glb->Jf_input_type == 7) {
				temporal_var = (pt_glb->time_device - pt_glb->dt) * (4.0 * PI * pt_glb->Jf_input_freq);
				temporal_var = temporal_var * temporal_var * temporal_var * (4.0 - temporal_var) * exp(-temporal_var);
				temporal_var = temporal_var / (dx * dy * dz);

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = temporal_var;
				}				
			}
		}
	}
	//	{
	//#pragma acc loop gang vector private(i,j,temporal_var, spatial_var)
	//		for (long int k = Jf_z1; k < Jf_z2; k++) {
	//			temporal_var = sin(2. * PI * pt_glb->Jf_test_freq * (pt_glb->time_device - pt_glb->dt));
	//			spatial_var = PI / Jf_thickness * (static_cast<double>(k - Jf_z1) + 0.5) * dz;
	//#pragma acc loop seq
	//			for (i = 0; i < nx; i++) {
	//#pragma acc loop seq
	//				for (j = 0; j < ny; j++) {
	//					if (pt_glb->material_cell(i, j, k) != 0) {
	//						if (pt_glb->Jf_test_type == 0) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//						}
	//						else if (pt_glb->Jf_test_type == 1) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//						}
	//						else if (pt_glb->Jf_test_type == 2) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
}

void EMdynamic_system::update_Jf_input_half() {
	unsigned long int i, j, k;
	double temporal_var, spatial_var;

#pragma acc parallel default(present) async(7)
	{
#pragma acc loop gang vector private(i,j,k, temporal_var, spatial_var)
		for (unsigned long id = 0; id < Jf_n; id++) {
			i = id / (Jf_ny * Jf_nz);
			j = (id - i * (Jf_ny * Jf_nz)) / Jf_nz;
			k = id - i * (Jf_ny * Jf_nz) - j * Jf_nz;

			i = i + pt_glb->Jfin_xi;
			j = j + pt_glb->Jfin_yi;
			k = k + pt_glb->Jfin_zi;

			if (pt_glb->Jf_input_type < 3) {
				temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
				spatial_var = PI / (static_cast<double>(Jf_nz) * dz) * (static_cast<double>(k - pt_glb->Jfin_zi) + 0.5) * dz;
				if (pt_glb->Jf_input_type == 0) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
				}
				else if (pt_glb->Jf_input_type == 1) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
				}
				else if (pt_glb->Jf_input_type == 2) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
				}
			}
			else if (pt_glb->Jf_input_type == 3) {
				if ((pt_glb->time_device - 0.5 * pt_glb->dt) <= static_cast<double>(pt_glb->num_Jf_cycles) / pt_glb->Jf_input_freq) {
					temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
				}
				else {
					temporal_var = 0.;
				}

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
			}
			else if (pt_glb->Jf_input_type == 4) {
				temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 5) {
				temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 6) {
				temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = cos(PI * pt_glb->Jf_rotate_xy / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * pt_glb->Jf_rotate_xy / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfx(i, j, k) = cos(PI * (pt_glb->Jf_rotate_xy + 90.) / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * (pt_glb->Jf_rotate_xy + 90.) / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
			}
			else if (pt_glb->Jf_input_type == 7) {
				temporal_var = (pt_glb->time_device - 0.5 * pt_glb->dt) * (4.0 * PI * pt_glb->Jf_input_freq);
				temporal_var = temporal_var * temporal_var * temporal_var * (4.0 - temporal_var) * exp(-temporal_var);
				temporal_var = temporal_var / (dx * dy * dz);

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = temporal_var;
				}				
			}
		}
	}

	//#pragma acc parallel default(present) async(7)
	//	{
	//#pragma acc loop gang vector private(i,j,temporal_var, spatial_var)
	//		for (long int k = Jf_z1; k < Jf_z2; k++) {
	//			temporal_var = sin(2. * PI * pt_glb->Jf_test_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
	//			spatial_var = PI / Jf_thickness * (static_cast<double>(k - Jf_z1) + 0.5) * dz;
	//#pragma acc loop seq
	//			for (i = 0; i < nx; i++) {
	//#pragma acc loop seq
	//				for (j = 0; j < ny; j++) {
	//					if (pt_glb->material_cell(i, j, k) != 0) {
	//						if (pt_glb->Jf_test_type == 0) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//						}
	//						else if (pt_glb->Jf_test_type == 1) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//						}
	//						else if (pt_glb->Jf_test_type == 2) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
}

void EMdynamic_system::update_Jf_input_full() {
	unsigned long int i, j, k;
	double temporal_var, spatial_var;

#pragma acc parallel default(present) async(7)
	{
#pragma acc loop gang vector private(i,j,k, temporal_var, spatial_var)
		for (unsigned long id = 0; id < Jf_n; id++) {
			i = id / (Jf_ny * Jf_nz);
			j = (id - i * (Jf_ny * Jf_nz)) / Jf_nz;
			k = id - i * (Jf_ny * Jf_nz) - j * Jf_nz;

			i = i + pt_glb->Jfin_xi;
			j = j + pt_glb->Jfin_yi;
			k = k + pt_glb->Jfin_zi;

			if (pt_glb->Jf_input_type < 3) {
				temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * pt_glb->time_device);
				spatial_var = PI / (static_cast<double>(Jf_nz) * dz) * (static_cast<double>(k - pt_glb->Jfin_zi) + 0.5) * dz;
				if (pt_glb->Jf_input_type == 0) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
				}
				else if (pt_glb->Jf_input_type == 1) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
				}
				else if (pt_glb->Jf_input_type == 2) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
				}
			}
			else if (pt_glb->Jf_input_type == 3) {
				if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_Jf_cycles) / pt_glb->Jf_input_freq) {
					temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device));
				}
				else {
					temporal_var = 0.;
				}

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
			}
			else if (pt_glb->Jf_input_type == 4) {
				temporal_var = pt_glb->time_device - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 5) {
				temporal_var = pt_glb->time_device - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 6) {
				temporal_var = pt_glb->time_device;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = cos(PI * pt_glb->Jf_rotate_xy / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * pt_glb->Jf_rotate_xy / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfx(i, j, k) = cos(PI * (pt_glb->Jf_rotate_xy + 90.) / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * (pt_glb->Jf_rotate_xy + 90.) / 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
			}
			else if (pt_glb->Jf_input_type == 7) {
				temporal_var = pt_glb->time_device * (4.0 * PI * pt_glb->Jf_input_freq);
				// printf("update_Jf_input_full: (time_device - dt) = %le, Jf_input_freq = %le, temporal_var = %le\n", pt_glb->time_device - pt_glb->dt, pt_glb->Jf_input_freq, temporal_var);

				temporal_var = temporal_var * temporal_var * temporal_var * (4.0 - temporal_var) * exp(-temporal_var);
				// printf("update_Jf_input_full: temporal_var = %le\n", temporal_var);

				temporal_var = temporal_var / (dx * dy * dz);
				// printf("update_Jf_input_full: Jf_input_amp = %le, dx = %le, dy = %le, dz = %le, temporal_var = %le\n", pt_glb->Jf_input_amp, dx, dy, dz, temporal_var);

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = temporal_var;
				}				
			}
		}
	}
	//	unsigned int i, j;
	//	double temporal_var, spatial_var;
	//
	//#pragma acc parallel default(present) async(7)
	//	{
	//#pragma acc loop gang vector private(i,j,temporal_var, spatial_var)
	//		for (long int k = Jf_z1; k < Jf_z2; k++) {
	//			temporal_var = sin(2. * PI * pt_glb->Jf_test_freq * pt_glb->time_device);
	//			spatial_var = PI / Jf_thickness * (static_cast<double>(k - Jf_z1) + 0.5) * dz;
	//#pragma acc loop seq
	//			for (i = 0; i < nx; i++) {
	//#pragma acc loop seq
	//				for (j = 0; j < ny; j++) {
	//					if (pt_glb->material_cell(i, j, k) != 0) {
	//						if (pt_glb->Jf_test_type == 0) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//						}
	//						else if (pt_glb->Jf_test_type == 1) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//						}
	//						else if (pt_glb->Jf_test_type == 2) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK1() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk1(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk1(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk1(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk1(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk1(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk1(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk1(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk1(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk1(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK2() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk2(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk2(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk2(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk2(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk2(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk2(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk2(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk2(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk2(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK3() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk3(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk3(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk3(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk3(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk3(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk3(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk3(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk3(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk3(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK4() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk4(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk4(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk4(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk4(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk4(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk4(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk4(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk4(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk4(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::update_Jp_RK1() {
#pragma acc loop gang vector
	for (long id = 0; id < n; id++) {
		Jpx_n1_store(id) = Jpx_n1(id) + dJpx_n1_rk1(id) * 0.5;
		Jpy_n1_store(id) = Jpy_n1(id) + dJpy_n1_rk1(id) * 0.5;
		Jpz_n1_store(id) = Jpz_n1(id) + dJpz_n1_rk1(id) * 0.5;
		Jpx_n2_store(id) = Jpx_n2(id) + dJpx_n2_rk1(id) * 0.5;
		Jpy_n2_store(id) = Jpy_n2(id) + dJpy_n2_rk1(id) * 0.5;
		Jpz_n2_store(id) = Jpz_n2(id) + dJpz_n2_rk1(id) * 0.5;
		Jpx_n3_store(id) = Jpx_n3(id) + dJpx_n3_rk1(id) * 0.5;
		Jpy_n3_store(id) = Jpy_n3(id) + dJpy_n3_rk1(id) * 0.5;
		Jpz_n3_store(id) = Jpz_n3(id) + dJpz_n3_rk1(id) * 0.5;
	}
}

#pragma acc routine gang
void EMdynamic_system::update_Jp_RK2() {
#pragma acc loop gang vector
	for (long id = 0; id < n; id++) {
		Jpx_n1_store(id) = Jpx_n1(id) + dJpx_n1_rk2(id) * 0.5;
		Jpy_n1_store(id) = Jpy_n1(id) + dJpy_n1_rk2(id) * 0.5;
		Jpz_n1_store(id) = Jpz_n1(id) + dJpz_n1_rk2(id) * 0.5;
		Jpx_n2_store(id) = Jpx_n2(id) + dJpx_n2_rk2(id) * 0.5;
		Jpy_n2_store(id) = Jpy_n2(id) + dJpy_n2_rk2(id) * 0.5;
		Jpz_n2_store(id) = Jpz_n2(id) + dJpz_n2_rk2(id) * 0.5;
		Jpx_n3_store(id) = Jpx_n3(id) + dJpx_n3_rk2(id) * 0.5;
		Jpy_n3_store(id) = Jpy_n3(id) + dJpy_n3_rk2(id) * 0.5;
		Jpz_n3_store(id) = Jpz_n3(id) + dJpz_n3_rk2(id) * 0.5;

		//printf("Jpx1: %le\t%le\t%ld\n", Jpx_n1(id), dJpx_n1_rk1(id), id);
		//printf("Jpy1: %le\t%le\t%ld\n", Jpy_n1(id), dJpy_n1_rk1(id), id);
		//printf("Jpz1: %le\t%le\t%ld\n", Jpz_n1(id), dJpz_n1_rk1(id), id);
		//printf("Jpx2: %le\t%le\t%ld\n", Jpx_n2(id), dJpx_n2_rk1(id), id);
		//printf("Jpy2: %le\t%le\t%ld\n", Jpy_n2(id), dJpy_n2_rk1(id), id);
		//printf("Jpz2: %le\t%le\t%ld\n", Jpz_n2(id), dJpz_n2_rk1(id), id);
		//printf("Jpx3: %le\t%le\t%ld\n", Jpx_n3(id), dJpx_n3_rk1(id), id);
		//printf("Jpy3: %le\t%le\t%ld\n", Jpy_n3(id), dJpy_n3_rk1(id), id);
		//printf("Jpz3: %le\t%le\t%ld\n", Jpz_n3(id), dJpz_n3_rk1(id), id);
	}
}

#pragma acc routine gang
void EMdynamic_system::update_Jp_RK3() {
#pragma acc loop gang vector
	for (long id = 0; id < n; id++) {
		Jpx_n1_store(id) = Jpx_n1(id) + dJpx_n1_rk3(id);
		Jpy_n1_store(id) = Jpy_n1(id) + dJpy_n1_rk3(id);
		Jpz_n1_store(id) = Jpz_n1(id) + dJpz_n1_rk3(id);
		Jpx_n2_store(id) = Jpx_n2(id) + dJpx_n2_rk3(id);
		Jpy_n2_store(id) = Jpy_n2(id) + dJpy_n2_rk3(id);
		Jpz_n2_store(id) = Jpz_n2(id) + dJpz_n2_rk3(id);
		Jpx_n3_store(id) = Jpx_n3(id) + dJpx_n3_rk3(id);
		Jpy_n3_store(id) = Jpy_n3(id) + dJpy_n3_rk3(id);
		Jpz_n3_store(id) = Jpz_n3(id) + dJpz_n3_rk3(id);
	}
}

#pragma acc routine gang
void EMdynamic_system::update_Jp() {
#pragma acc loop gang vector
	for (long int id = 0; id < n; id++) {
		Jpx_n1(id) = Jpx_n1(id) + dJpx_n1_rk1(id) / 6. + dJpx_n1_rk2(id) / 3. + dJpx_n1_rk3(id) / 3. + dJpx_n1_rk4(id) / 6.;
		Jpx_n1_store(id) = Jpx_n1(id);
		Jpy_n1(id) = Jpy_n1(id) + dJpy_n1_rk1(id) / 6. + dJpy_n1_rk2(id) / 3. + dJpy_n1_rk3(id) / 3. + dJpy_n1_rk4(id) / 6.;
		Jpy_n1_store(id) = Jpy_n1(id);
		Jpz_n1(id) = Jpz_n1(id) + dJpz_n1_rk1(id) / 6. + dJpz_n1_rk2(id) / 3. + dJpz_n1_rk3(id) / 3. + dJpz_n1_rk4(id) / 6.;
		Jpz_n1_store(id) = Jpz_n1(id);
		Jpx_n2(id) = Jpx_n2(id) + dJpx_n2_rk1(id) / 6. + dJpx_n2_rk2(id) / 3. + dJpx_n2_rk3(id) / 3. + dJpx_n2_rk4(id) / 6.;
		Jpx_n2_store(id) = Jpx_n2(id);
		Jpy_n2(id) = Jpy_n2(id) + dJpy_n2_rk1(id) / 6. + dJpy_n2_rk2(id) / 3. + dJpy_n2_rk3(id) / 3. + dJpy_n2_rk4(id) / 6.;
		Jpy_n2_store(id) = Jpy_n2(id);
		Jpz_n2(id) = Jpz_n2(id) + dJpz_n2_rk1(id) / 6. + dJpz_n2_rk2(id) / 3. + dJpz_n2_rk3(id) / 3. + dJpz_n2_rk4(id) / 6.;
		Jpz_n2_store(id) = Jpz_n2(id);
		Jpx_n3(id) = Jpx_n3(id) + dJpx_n3_rk1(id) / 6. + dJpx_n3_rk2(id) / 3. + dJpx_n3_rk3(id) / 3. + dJpx_n3_rk4(id) / 6.;
		Jpx_n3_store(id) = Jpx_n3(id);
		Jpy_n3(id) = Jpy_n3(id) + dJpy_n3_rk1(id) / 6. + dJpy_n3_rk2(id) / 3. + dJpy_n3_rk3(id) / 3. + dJpy_n3_rk4(id) / 6.;
		Jpy_n3_store(id) = Jpy_n3(id);
		Jpz_n3(id) = Jpz_n3(id) + dJpz_n3_rk1(id) / 6. + dJpz_n3_rk2(id) / 3. + dJpz_n3_rk3(id) / 3. + dJpz_n3_rk4(id) / 6.;
		Jpz_n3_store(id) = Jpz_n3(id);
	}
}
