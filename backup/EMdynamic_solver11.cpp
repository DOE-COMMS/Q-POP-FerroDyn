#include "EMdynamic_system.h"

//void EMdynamic_system::get_dE() {
//	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
//	double dP;
//
//	long int i, j, k;
//	unsigned int mat_type;
//	material* mat;
//	long int idx1, idy1, idz1;
//	long int idx2, idy2, idz2;
//	long int idx3, idy3, idz3;
//	long int idx4, idy4, idz4;
//
//	double er;
//	double cond;
//	double Jf, Jp, Jishe;
//
//	if (pt_glb->if_Jf_test == true) {
//		update_Jf_test();
//	}
//#pragma acc wait
//
//	//-----------X component-----------//
//#pragma acc parallel loop gang vector \
//private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
//er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k) default(present) async(1)
//	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
//		i = id / ((ny + 1) * (nz + 1));
//		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
//		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);
//
//		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
//		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
//			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
//		}
//		else {
//			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
//		}
//		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
//			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
//		}
//		else {
//			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
//		}
//
//		dHzdy = (DHz_em(idx1, idy3, idz4) - DHz_em(idx1, idy1, idz4)) / dy;
//		dHydz = (DHy_em(idx1, idy3, idz2) - DHy_em(idx1, idy3, idz1)) / dz;
//
//		Jf = 0.; Jp = 0.; Jishe = 0.;
//		dP = 0.; er = 0.; cond = 0.;
//		//----------Polarization 1------//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP + pt_fe->dpx_glb(idx1, idy1, idz1);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx1, idy1, idz1);
//		Jp = Jp + Jpx_n1(idx1, idy1, idz1) + Jpx_n2(idx1, idy1, idz1) + Jpx_n3(idx1, idy1, idz1);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);
//
//		//----------Polarization 2------//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpx_glb(idx2, idy2, idz2);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx2, idy2, idz2);
//		Jp = Jp + Jpx_n1(idx2, idy2, idz2) + Jpx_n2(idx2, idy2, idz2) + Jpx_n3(idx2, idy2, idz2);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);
//
//		//----------Polarization 3------//
//		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpx_glb(idx3, idy3, idz3);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx3, idy3, idz3);
//		Jp = Jp + Jpx_n1(idx3, idy3, idz3) + Jpx_n2(idx3, idy3, idz3) + Jpx_n3(idx3, idy3, idz3);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);
//
//		//----------Polarization 4------//
//		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpx_glb(idx4, idy4, idz4);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx4, idy4, idz4);
//		Jp = Jp + Jpx_n1(idx4, idy4, idz4) + Jpx_n2(idx4, idy4, idz4) + Jpx_n3(idx4, idy4, idz4);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);
//
//		Jf = Jf / 4.;
//		Jp = Jp / 4.;
//		Jishe = Jishe / 4.;
//		dP = dP / 4.; er = er / 4.; cond = cond / 4.;
//
//		dDEx_em(i, j, k) = pt_glb->dt / e0 / er * \
//			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em(i, j, k) - dP/ pt_glb->dt);
//	}
//	//-----------END X component-----------//
//
//	//-----------Y component-----------//
//#pragma acc parallel loop gang vector \
//private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
//Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k) default(present) async(2)
//	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
//		i = id / ((ny) * (nz + 1));
//		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
//		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);
//
//		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
//			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
//		}
//		else {
//			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
//		}
//
//		idy1 = j; idy2 = j; idy3 = j; idy4 = j;
//
//		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
//			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
//		}
//		else {
//			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
//		}
//
//		dHxdz = (DHx_em(idx3, idy1, idz4) - DHx_em(idx3, idy1, idz3)) / dz;
//		dHzdx = (DHz_em(idx3, idy1, idz2) - DHz_em(idx1, idy1, idz2)) / dx;
//
//		Jf = 0.; Jp = 0.; Jishe = 0.;
//		dP = 0.; er = 0.; cond = 0.;
//		//----------Polarization 1------//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP =dP+ pt_fe->dpy_glb(idx1, idy1, idz1);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx1, idy1, idz1);
//		Jp = Jp + Jpy_n1(idx1, idy1, idz1) + Jpy_n2(idx1, idy1, idz1) + Jpy_n3(idx1, idy1, idz1);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);
//
//		//----------Polarization 2------//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpy_glb(idx2, idy2, idz2);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx2, idy2, idz2);
//		Jp = Jp + Jpy_n1(idx2, idy2, idz2) + Jpy_n2(idx2, idy2, idz2) + Jpy_n3(idx2, idy2, idz2);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);
//
//		//----------Polarization 3------//
//		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpy_glb(idx3, idy3, idz3);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx3, idy3, idz3);
//		Jp = Jp + Jpy_n1(idx3, idy3, idz3) + Jpy_n2(idx3, idy3, idz3) + Jpy_n3(idx3, idy3, idz3);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);
//
//		//----------Polarization 4------//
//		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpy_glb(idx4, idy4, idz4);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx4, idy4, idz4);
//		Jp = Jp + Jpy_n1(idx4, idy4, idz4) + Jpy_n2(idx4, idy4, idz4) + Jpy_n3(idx4, idy4, idz4);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);
//
//		Jf = Jf / 4.;
//		Jp = Jp / 4.;
//		Jishe = Jishe / 4.;
//		dP = dP / 4.; er = er / 4.; cond = cond / 4.;
//
//		dDEy_em(i, j, k) = pt_glb->dt / e0 / er * \
//			(dHxdz - dHzdx - Jf - Jp - Jishe -cond * DEy_em(i, j, k) - dP/ pt_glb->dt);
//	}
//
//	//-----------END Y component-----------//
//
//	//-----------Z component-----------//
//#pragma acc parallel loop gang vector \
//private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
//er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k) default(present) async(3)
//	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
//		i = id / ((ny + 1) * (nz));
//		j = (id - i * ((ny + 1) * (nz))) / (nz);
//		k = id - i * ((ny + 1) * (nz)) - j * (nz);
//
//		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
//			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
//		}
//		else {
//			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
//		}
//
//		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
//			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
//		}
//		else {
//			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
//		}
//
//		idz1 = k; idz2 = k; idz3 = k; idz4 = k;
//
//		dHydx = (DHy_em(idx3, idy2, idz1) - DHy_em(idx1, idy2, idz1)) / dx;
//		dHxdy = (DHx_em(idx3, idy2, idz1) - DHx_em(idx3, idy1, idz1)) / dy;
//
//		Jf = 0.; Jp = 0.;
//		dP = 0.; er = 0.; cond = 0.;
//		//----------Polarization 1------//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+pt_fe->dpz_glb(idx1, idy1, idz1);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx1, idy1, idz1);
//		Jp = Jp + Jpz_n1(idx1, idy1, idz1) + Jpz_n2(idx1, idy1, idz1) + Jpz_n3(idx1, idy1, idz1);
//
//		//----------Polarization 2------//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpz_glb(idx2, idy2, idz2);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx2, idy2, idz2);
//		Jp = Jp + Jpz_n1(idx2, idy2, idz2) + Jpz_n2(idx2, idy2, idz2) + Jpz_n3(idx2, idy2, idz2);
//
//		//----------Polarization 3------//
//		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP =dP+ pt_fe->dpz_glb(idx3, idy3, idz3);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx3, idy3, idz3);
//		Jp = Jp + Jpz_n1(idx3, idy3, idz3) + Jpz_n2(idx3, idy3, idz3) + Jpz_n3(idx3, idy3, idz3);
//
//		//----------Polarization 4------//
//		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP =dP+ pt_fe->dpz_glb(idx4, idy4, idz4);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx4, idy4, idz4);
//		Jp = Jp + Jpz_n1(idx4, idy4, idz4) + Jpz_n2(idx4, idy4, idz4) + Jpz_n3(idx4, idy4, idz4);
//
//		Jf = Jf / 4.;
//		Jp = Jp / 4.;
//		dP = dP / 4.; er = er / 4.; cond = cond / 4.;
//
//		dDEz_em(i, j, k) = pt_glb->dt / e0 / er * \
//			(dHydx - dHxdy - Jf - Jp - cond * DEz_em(i, j, k) -dP/ pt_glb->dt);
//	}
//	//-----------END Z component-----------//
//}

void EMdynamic_system::get_dE_RK1() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int xidx1, xidy1, xidz1, yidx1, yidy1, yidz1, zidx1, zidy1, zidz1;
	long int xidx2, xidy2, xidz2, yidx2, yidy2, yidz2, zidx2, zidy2, zidz2;
	long int xidx3, xidy3, xidz3, yidx3, yidy3, yidz3, zidx3, zidy3, zidz3;
	long int xidx4, xidy4, xidz4, yidx4, yidy4, yidz4, zidx4, zidy4, zidz4;

	double xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33;
	double yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33;
	double zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33;

	double xDenominator, yDenominator, zDenominator;

	double xinverse_er11, xinverse_er12, xinverse_er13;
	double xinverse_er21, xinverse_er22, xinverse_er23;
	double xinverse_er31, xinverse_er32, xinverse_er33;
	double yinverse_er11, yinverse_er12, yinverse_er13;  
	double yinverse_er21, yinverse_er22, yinverse_er23;
	double yinverse_er31, yinverse_er32, yinverse_er33;
	double zinverse_er11, zinverse_er12, zinverse_er13;
	double zinverse_er21, zinverse_er22, zinverse_er23;
	double zinverse_er31, zinverse_er32, zinverse_er33;



	double xcond, ycond, zcond;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double Px_count, Py_count, Pz_count;

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
private(xidx1, xidy1, xidz1, xidx2, xidy2, xidz2, xidx3, xidy3, xidz3, xidx4, xidy4, xidz4,\
yidx1, yidy1, yidz1, yidx2, yidy2, yidz2, yidx3, yidy3, yidz3, yidx4, yidy4, yidz4,\
zidx1, zidy1, zidz1, zidx2, zidy2, zidz2, zidx3, zidy3, zidz3, zidx4, zidy4, zidz4,\
xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33,\
yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33,\
zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33,\
xDenominator,\
xinverse_er11, xinverse_er12, xinverse_er13,\
xinverse_er21, xinverse_er22, xinverse_er23,\
xinverse_er31, xinverse_er32, xinverse_er33,\
yDenominator,\
yinverse_er11, yinverse_er12, yinverse_er13,\
yinverse_er21, yinverse_er22, yinverse_er23,\
yinverse_er31, yinverse_er32, yinverse_er33,\
zDenominator,\
zinverse_er11, zinverse_er12, zinverse_er13,\
zinverse_er21, zinverse_er22, zinverse_er23,\
zinverse_er31, zinverse_er32, zinverse_er33,\
xcond, ycond, zcond,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,\
Px_count,Py_count,Pz_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);
		
		
		// Take the x component
		xidx1 = i; xidx2 = i; xidx3 = i; xidx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			xidy1 = ny - 1; xidy2 = ny - 1; xidy3 = 0; xidy4 = 0;
		}
		else {
			xidy1 = j - 1; xidy2 = j - 1; xidy3 = j; xidy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			xidz1 = nz - 1; xidz2 = 0; xidz3 = nz - 1; xidz4 = 0;
		}
		else {
			xidz1 = k - 1; xidz2 = k; xidz3 = k - 1; xidz4 = k;
		}

		dHzdy = (DHz_em_store(xidx1, xidy3, xidz4) - DHz_em_store(xidx1, xidy1, xidz4)) / dy;
		dHydz = (DHy_em_store(xidx1, xidy3, xidz2) - DHy_em_store(xidx1, xidy3, xidz1)) / dz;

		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.; 

		
		// Take the y component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			yidx1 = nx - 1; yidx2 = nx - 1; yidx3 = 0; yidx4 = 0;
		}
		else {
			yidx1 = i - 1; yidx2 = i - 1; yidx3 = i; yidx4 = i;
		}

		yidy1 = j; yidy2 = j; yidy3 = j; yidy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			yidz1 = nz - 1; yidz2 = 0; yidz3 = nz - 1; yidz4 = 0;
		}
		else {
			yidz1 = k - 1; yidz2 = k; yidz3 = k - 1; yidz4 = k;
		}

		dHxdz = (DHx_em_store(yidx3, yidy1, yidz4) - DHx_em_store(yidx3, yidy1, yidz3)) / dz;
		dHzdx = (DHz_em_store(yidx3, yidy1, yidz2) - DHz_em_store(yidx1, yidy1, yidz2)) / dx;

		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.; 

		
		// Take the z component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			zidx1 = nx - 1; zidx2 = nx - 1; zidx3 = 0; zidx4 = 0;
		}
		else {
			zidx1 = i - 1; zidx2 = i - 1; zidx3 = i; zidx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			zidy1 = ny - 1; zidy2 = 0; zidy3 = ny - 1; zidy4 = 0;
		}
		else {
			zidy1 = j - 1; zidy2 = j; zidy3 = j - 1; zidy4 = j;
		}

		zidz1 = k; zidz2 = k; zidz3 = k; zidz4 = k;

		dHydx = (DHy_em_store(zidx3, zidy2, zidz1) - DHy_em_store(zidx1, zidy2, zidz1)) / dx;
		dHxdy = (DHx_em_store(zidx3, zidy2, zidz1) - DHx_em_store(zidx3, zidy1, zidz1)) / dy;

		Jfz = 0.; Jpz = 0.;
		dPz = 0.;
		
		
		
		xcond = 0.;
		ycond = 0.;
		zcond = 0.;

		xer11 = 0.; xer12 = 0.; xer13 = 0.;
		xer21 = 0.; xer22 = 0.; xer23 = 0.;
		xer31 = 0.; xer32 = 0.; xer33 = 0.;

		yer11 = 0.; yer12 = 0.; yer13 = 0.;
		yer21 = 0.; yer22 = 0.; yer23 = 0.;
		yer31 = 0.; yer32 = 0.; yer33 = 0.;

		zer11 = 0.; zer12 = 0.; zer13 = 0.;
		zer21 = 0.; zer22 = 0.; zer23 = 0.;
		zer31 = 0.; zer32 = 0.; zer33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		Px_count = 0.; Py_count = 0.; Pz_count = 0.;



		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(xidx1, xidy1, xidz1);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(xidx1, xidy1, xidz1); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx1, xidy1, xidz1);
		Jpx = Jpx + Jpx_n1_store(xidx1, xidy1, xidz1) + Jpx_n2_store(xidx1, xidy1, xidz1) + Jpx_n3_store(xidx1, xidy1, xidz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx1, xidy1, xidz1);


		mat_type = pt_glb->material_cell(yidx1, yidy1, yidz1);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk1(yidx1, yidy1, yidz1); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx1, yidy1, yidz1);
		Jpy = Jpy + Jpy_n1_store(yidx1, yidy1, yidz1) + \
			Jpy_n2_store(yidx1, yidy1, yidz1) + Jpy_n3_store(yidx1, yidy1, yidz1);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx1, yidy1, yidz1);


		mat_type = pt_glb->material_cell(zidx1, zidy1, zidz1);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk1(zidx1, zidy1, zidz1); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx1, zidy1, zidz1);
		Jpz = Jpz + Jpz_n1_store(zidx1, zidy1, zidz1) + \
			Jpz_n2_store(zidx1, zidy1, zidz1) + Jpz_n3_store(zidx1, zidy1, zidz1);




		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(xidx2, xidy2, xidz2);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(xidx2, xidy2, xidz2); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx2, xidy2, xidz2);
		Jpx = Jpx + Jpx_n1_store(xidx2, xidy2, xidz2) + Jpx_n2_store(xidx2, xidy2, xidz2) + Jpx_n3_store(xidx2, xidy2, xidz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx2, xidy2, xidz2);


		mat_type = pt_glb->material_cell(yidx2, yidy2, yidz2);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk1(yidx2, yidy2, yidz2); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx2, yidy2, yidz2);
		Jpy = Jpy + Jpy_n1_store(yidx2, yidy2, yidz2) + \
			Jpy_n2_store(yidx2, yidy2, yidz2) + Jpy_n3_store(yidx2, yidy2, yidz2);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx2, yidy2, yidz2);


		mat_type = pt_glb->material_cell(zidx2, zidy2, zidz2);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk1(zidx2, zidy2, zidz2); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx2, zidy2, zidz2);
		Jpz = Jpz + Jpz_n1_store(zidx2, zidy2, zidz2) + \
			Jpz_n2_store(zidx2, zidy2, zidz2) + Jpz_n3_store(zidx2, zidy2, zidz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(xidx3, xidy3, xidz3);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(xidx3, xidy3, xidz3); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx3, xidy3, xidz3);
		Jpx = Jpx + Jpx_n1_store(xidx3, xidy3, xidz3) + \
			Jpx_n2_store(xidx3, xidy3, xidz3) + Jpx_n3_store(xidx3, xidy3, xidz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx3, xidy3, xidz3);


		mat_type = pt_glb->material_cell(yidx3, yidy3, yidz3);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk1(yidx3, yidy3, yidz3); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx3, yidy3, yidz3);
		Jpy = Jpy + Jpy_n1_store(yidx3, yidy3, yidz3) + \
			Jpy_n2_store(yidx3, yidy3, yidz3) + Jpy_n3_store(yidx3, yidy3, yidz3);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx3, yidy3, yidz3);


		mat_type = pt_glb->material_cell(zidx3, zidy3, zidz3);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk1(zidx3, zidy3, zidz3); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx3, zidy3, zidz3);
		Jpz = Jpz + Jpz_n1_store(zidx3, zidy3, zidz3) + \
			Jpz_n2_store(zidx3, zidy3, zidz3) + Jpz_n3_store(zidx3, zidy3, zidz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(xidx4, xidy4, xidz4);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk1(xidx4, xidy4, xidz4); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx4, xidy4, xidz4);
		Jpx = Jpx + Jpx_n1_store(xidx4, xidy4, xidz4) + \
			Jpx_n2_store(xidx4, xidy4, xidz4) + Jpx_n3_store(xidx4, xidy4, xidz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx4, xidy4, xidz4);


		mat_type = pt_glb->material_cell(yidx4, yidy4, yidz4);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk1(yidx4, yidy4, yidz4); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx4, yidy4, yidz4);
		Jpy = Jpy + Jpy_n1_store(yidx4, yidy4, yidz4) + \
			Jpy_n2_store(yidx4, yidy4, yidz4) + Jpy_n3_store(yidx4, yidy4, yidz4);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx4, yidy4, yidz4);


		mat_type = pt_glb->material_cell(zidx4, zidy4, zidz4);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk1(zidx4, zidy4, zidz4); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx4, zidy4, zidz4);
		Jpz = Jpz + Jpz_n1_store(zidx4, zidy4, zidz4) + \
			Jpz_n2_store(zidx4, zidy4, zidz4) + Jpz_n3_store(zidx4, zidy4, zidz4);


		if (Px_count < 0.5) {
			Px_count = 1.;
		}
		if (Py_count < 0.5) {
			Py_count = 1.;
		}
		if (Pz_count < 0.5) {
			Pz_count = 1.;
		}

		Jfx = Jfx / 4.;	Jfy = Jfy / 4.;	Jfz = Jfz / 4.;
		Jpx = Jpx / 4.;	Jpy = Jpy / 4.;	Jpz = Jpz / 4.;	
		Jishex = Jishex / 4.;	Jishey = Jishey / 4.;	
		dPx = dPx / Px_count;	dPy = dPy / Py_count;	dPz = dPz / Pz_count;
		
		xcond = xcond / 4.;	ycond = ycond / 4.;	zcond = zcond / 4.;

		xer11 = xer11 / 4.; xer12 = xer12 / 4.; xer13 = xer13 / 4.;
		xer21 = xer21 / 4.; xer22 = xer22 / 4.; xer23 = xer23 / 4.;
		xer31 = xer31 / 4.; xer32 = xer32 / 4.; xer33 = xer33 / 4.;

		yer11 = yer11 / 4.; yer12 = yer12 / 4.; yer13 = yer13 / 4.;
		yer21 = yer21 / 4.; yer22 = yer22 / 4.; yer23 = yer23 / 4.;
		yer31 = yer31 / 4.; yer32 = yer32 / 4.; yer33 = yer33 / 4.;

		zer11 = zer11 / 4.; zer12 = zer12 / 4.; zer13 = zer13 / 4.;
		zer21 = zer21 / 4.; zer22 = zer22 / 4.; zer23 = zer23 / 4.;
		zer31 = zer31 / 4.; zer32 = zer32 / 4.; zer33 = zer33 / 4.;

		xDenominator = -xer13*xer22*xer31 + xer12*xer23*xer31 + xer13*xer21*xer32 -\
			xer11*xer23*xer32 - xer12*xer21*xer33 + xer11*xer22*xer33;
		yDenominator = -yer13*yer22*yer31 + yer12*yer23*yer31 + yer13*yer21*yer32 -\
			yer11*yer23*yer32 - yer12*yer21*yer33 + yer11*yer22*yer33;
		zDenominator = -zer13*zer22*zer31 + zer12*zer23*zer31 + zer13*zer21*zer32 -\
			zer11*zer23*zer32 - zer12*zer21*zer33 + zer11*zer22*zer33;

		xinverse_er11 = (-xer23*xer32 + xer22*xer33) / xDenominator;
		xinverse_er12 = ( xer13*xer32 - xer12*xer33) / xDenominator;
		xinverse_er13 = (-xer13*xer22 + xer12*xer23) / xDenominator;
		xinverse_er21 = ( xer23*xer31 - xer21*xer33) / xDenominator;
		xinverse_er22 = (-xer13*xer31 + xer11*xer33) / xDenominator;
		xinverse_er23 = ( xer13*xer21 - xer11*xer23) / xDenominator;
		xinverse_er31 = (-xer22*xer31 + xer21*xer32) / xDenominator;
		xinverse_er32 = ( xer12*xer31 - xer11*xer32) / xDenominator;
		xinverse_er33 = (-xer12*xer21 + xer11*xer22) / xDenominator;


		yinverse_er11 = (-yer23*yer32 + yer22*yer33) / yDenominator;
		yinverse_er12 = ( yer13*yer32 - yer12*yer33) / yDenominator;
		yinverse_er13 = (-yer13*yer22 + yer12*yer23) / yDenominator;
		yinverse_er21 = ( yer23*yer31 - yer21*yer33) / yDenominator;
		yinverse_er22 = (-yer13*yer31 + yer11*yer33) / yDenominator;
		yinverse_er23 = ( yer13*yer21 - yer11*yer23) / yDenominator;
		yinverse_er31 = (-yer22*yer31 + yer21*yer32) / yDenominator;
		yinverse_er32 = ( yer12*yer31 - yer11*yer32) / yDenominator;
		yinverse_er33 = (-yer12*yer21 + yer11*yer22) / yDenominator;


		zinverse_er11 = (-zer23*zer32 + zer22*zer33) / zDenominator;
		zinverse_er12 = ( zer13*zer32 - zer12*zer33) / zDenominator;
		zinverse_er13 = (-zer13*zer22 + zer12*zer23) / zDenominator;
		zinverse_er21 = ( zer23*zer31 - zer21*zer33) / zDenominator;
		zinverse_er22 = (-zer13*zer31 + zer11*zer33) / zDenominator;
		zinverse_er23 = ( zer13*zer21 - zer11*zer23) / zDenominator;
		zinverse_er31 = (-zer22*zer31 + zer21*zer32) / zDenominator;
		zinverse_er32 = ( zer12*zer31 - zer11*zer32) / zDenominator;
		zinverse_er33 = (-zer12*zer21 + zer11*zer22) / zDenominator;



		dDEx_em_rk1(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er11 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er12 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er13 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEy_em_rk1(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er21 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er22 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er23 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEz_em_rk1(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er31 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er32 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er33 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK
	}
	//-----------END X Y Z component-----------//









// 	//-----------X component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k, P_count) default(present) async(1)
// 	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
// 		i = id / ((ny + 1) * (nz + 1));
// 		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

// 		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
// 		}
// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
// 		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; cond = 0.;
// 		er11 = 0.; er12 = 0.; er13 = 0.;
// 		er21 = 0.; er22 = 0.; er23 = 0.;
// 		er31 = 0.; er32 = 0.; er33 = 0.;

// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er11 = er11 + 1.; 	er12 = er12 + 0.; 	er13 = er13 + 0.;
// 			er21 = er21 + 0.; 	er22 = er22 + 1.; 	er23 = er23 + 0.;
// 			er31 = er31 + 0.; 	er32 = er32 + 0.; 	er33 = er33 + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk1(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			//if (mat->if_spin_pump == true) {
// 			//	Ji_count = Ji_count + 1.;
// 			//}
// 			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
// 			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
// 			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx1, idy1, idz1);
// 		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er11 = er11 + 1.; 	er12 = er12 + 0.; 	er13 = er13 + 0.;
// 			er21 = er21 + 0.; 	er22 = er22 + 1.; 	er23 = er23 + 0.;
// 			er31 = er31 + 0.; 	er32 = er32 + 0.; 	er33 = er33 + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk1(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
// 			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
// 			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx2, idy2, idz2);
// 		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er11 = er11 + 1.; 	er12 = er12 + 0.; 	er13 = er13 + 0.;
// 			er21 = er21 + 0.; 	er22 = er22 + 1.; 	er23 = er23 + 0.;
// 			er31 = er31 + 0.; 	er32 = er32 + 0.; 	er33 = er33 + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk1(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
// 			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
// 			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx3, idy3, idz3);
// 		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er11 = er11 + 1.; 	er12 = er12 + 0.; 	er13 = er13 + 0.;
// 			er21 = er21 + 0.; 	er22 = er22 + 1.; 	er23 = er23 + 0.;
// 			er31 = er31 + 0.; 	er32 = er32 + 0.; 	er33 = er33 + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk1(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er11 = er11 + mat->r_permittivity11; er12 = er12 + mat->r_permittivity12; er13 = er13 + mat->r_permittivity13;
// 			er21 = er21 + mat->r_permittivity21; er22 = er22 + mat->r_permittivity22; er23 = er23 + mat->r_permittivity23;
// 			er31 = er31 + mat->r_permittivity31; er32 = er32 + mat->r_permittivity32; er33 = er33 + mat->r_permittivity33;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx4, idy4, idz4);
// 		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; cond = cond / 4.;

// 		er11 = er11 / 4.; er12 = er12 / 4.; er13 = er13 / 4.;
// 		er21 = er21 / 4.; er22 = er22 / 4.; er23 = er23 / 4.;
// 		er31 = er31 / 4.; er32 = er32 / 4.; er33 = er33 / 4.;

// 		dDEx_em_rk1(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END X component-----------//

// 	//-----------Y component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
// 	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
// 		i = id / ((ny) * (nz + 1));
// 		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
// 		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk1(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx1, idy1, idz1);
// 		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk1(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx2, idy2, idz2);
// 		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk1(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx3, idy3, idz3);
// 		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk1(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx4, idy4, idz4);
// 		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEy_em_rk1(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}

// 	//-----------END Y component-----------//

// 	//-----------Z component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
// 	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
// 		i = id / ((ny + 1) * (nz));
// 		j = (id - i * ((ny + 1) * (nz))) / (nz);
// 		k = id - i * ((ny + 1) * (nz)) - j * (nz);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
// 		}

// 		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

// 		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
// 		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

// 		Jf = 0.; Jp = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk1(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx1, idy1, idz1);
// 		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk1(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx2, idy2, idz2);
// 		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk1(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx3, idy3, idz3);
// 		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk1(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx4, idy4, idz4);
// 		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEz_em_rk1(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END Z component-----------//
}

void EMdynamic_system::get_dE_RK2() {
	// double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	// double dP;

	// long int i, j, k;
	// unsigned int mat_type;
	// material* mat;
	// long int idx1, idy1, idz1;
	// long int idx2, idy2, idz2;
	// long int idx3, idy3, idz3;
	// long int idx4, idy4, idz4;

	// double er;
	// double cond;
	// double Jf, Jp, Jishe;
	// //double Jf_count, Jp_count, Ji_count;
	// double P_count;

	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int xidx1, xidy1, xidz1, yidx1, yidy1, yidz1, zidx1, zidy1, zidz1;
	long int xidx2, xidy2, xidz2, yidx2, yidy2, yidz2, zidx2, zidy2, zidz2;
	long int xidx3, xidy3, xidz3, yidx3, yidy3, yidz3, zidx3, zidy3, zidz3;
	long int xidx4, xidy4, xidz4, yidx4, yidy4, yidz4, zidx4, zidy4, zidz4;

	double xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33;
	double yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33;
	double zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33;

	double xDenominator, yDenominator, zDenominator;

	double xinverse_er11, xinverse_er12, xinverse_er13;
	double xinverse_er21, xinverse_er22, xinverse_er23;
	double xinverse_er31, xinverse_er32, xinverse_er33;
	double yinverse_er11, yinverse_er12, yinverse_er13;  
	double yinverse_er21, yinverse_er22, yinverse_er23;
	double yinverse_er31, yinverse_er32, yinverse_er33;
	double zinverse_er11, zinverse_er12, zinverse_er13;
	double zinverse_er21, zinverse_er22, zinverse_er23;
	double zinverse_er31, zinverse_er32, zinverse_er33;



	double xcond, ycond, zcond;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double Px_count, Py_count, Pz_count;





	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_half(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK2();
	}

#pragma acc wait
// 	//-----------X component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k,P_count) default(present) async(1)
// 	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
// 		i = id / ((ny + 1) * (nz + 1));
// 		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

// 		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
// 		}
// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
// 		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk2(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx1, idy1, idz1);
// 		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk2(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx2, idy2, idz2);
// 		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk2(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx3, idy3, idz3);
// 		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk2(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx4, idy4, idz4);
// 		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEx_em_rk2(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END X component-----------//

// 	//-----------Y component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
// 	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
// 		i = id / ((ny) * (nz + 1));
// 		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
// 		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk2(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx1, idy1, idz1);
// 		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk2(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx2, idy2, idz2);
// 		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk2(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx3, idy3, idz3);
// 		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk2(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx4, idy4, idz4);
// 		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEy_em_rk2(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}

// 	//-----------END Y component-----------//

// 	//-----------Z component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
// 	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
// 		i = id / ((ny + 1) * (nz));
// 		j = (id - i * ((ny + 1) * (nz))) / (nz);
// 		k = id - i * ((ny + 1) * (nz)) - j * (nz);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
// 		}

// 		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

// 		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
// 		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

// 		Jf = 0.; Jp = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk2(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx1, idy1, idz1);
// 		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk2(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx2, idy2, idz2);
// 		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk2(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx3, idy3, idz3);
// 		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk2(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx4, idy4, idz4);
// 		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEz_em_rk2(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END Z component-----------//


	//-----------X Y Z component-----------//
#pragma acc parallel loop gang vector \
private(xidx1, xidy1, xidz1, xidx2, xidy2, xidz2, xidx3, xidy3, xidz3, xidx4, xidy4, xidz4,\
yidx1, yidy1, yidz1, yidx2, yidy2, yidz2, yidx3, yidy3, yidz3, yidx4, yidy4, yidz4,\
zidx1, zidy1, zidz1, zidx2, zidy2, zidz2, zidx3, zidy3, zidz3, zidx4, zidy4, zidz4,\
xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33,\
yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33,\
zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33,\
xDenominator,\
xinverse_er11, xinverse_er12, xinverse_er13,\
xinverse_er21, xinverse_er22, xinverse_er23,\
xinverse_er31, xinverse_er32, xinverse_er33,\
yDenominator,\
yinverse_er11, yinverse_er12, yinverse_er13,\
yinverse_er21, yinverse_er22, yinverse_er23,\
yinverse_er31, yinverse_er32, yinverse_er33,\
zDenominator,\
zinverse_er11, zinverse_er12, zinverse_er13,\
zinverse_er21, zinverse_er22, zinverse_er23,\
zinverse_er31, zinverse_er32, zinverse_er33,\
xcond, ycond, zcond,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,\
Px_count,Py_count,Pz_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);
		
		
		// Take the x component
		xidx1 = i; xidx2 = i; xidx3 = i; xidx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			xidy1 = ny - 1; xidy2 = ny - 1; xidy3 = 0; xidy4 = 0;
		}
		else {
			xidy1 = j - 1; xidy2 = j - 1; xidy3 = j; xidy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			xidz1 = nz - 1; xidz2 = 0; xidz3 = nz - 1; xidz4 = 0;
		}
		else {
			xidz1 = k - 1; xidz2 = k; xidz3 = k - 1; xidz4 = k;
		}

		dHzdy = (DHz_em_store(xidx1, xidy3, xidz4) - DHz_em_store(xidx1, xidy1, xidz4)) / dy;
		dHydz = (DHy_em_store(xidx1, xidy3, xidz2) - DHy_em_store(xidx1, xidy3, xidz1)) / dz;

		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.; 

		
		// Take the y component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			yidx1 = nx - 1; yidx2 = nx - 1; yidx3 = 0; yidx4 = 0;
		}
		else {
			yidx1 = i - 1; yidx2 = i - 1; yidx3 = i; yidx4 = i;
		}

		yidy1 = j; yidy2 = j; yidy3 = j; yidy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			yidz1 = nz - 1; yidz2 = 0; yidz3 = nz - 1; yidz4 = 0;
		}
		else {
			yidz1 = k - 1; yidz2 = k; yidz3 = k - 1; yidz4 = k;
		}

		dHxdz = (DHx_em_store(yidx3, yidy1, yidz4) - DHx_em_store(yidx3, yidy1, yidz3)) / dz;
		dHzdx = (DHz_em_store(yidx3, yidy1, yidz2) - DHz_em_store(yidx1, yidy1, yidz2)) / dx;

		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.; 

		
		// Take the z component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			zidx1 = nx - 1; zidx2 = nx - 1; zidx3 = 0; zidx4 = 0;
		}
		else {
			zidx1 = i - 1; zidx2 = i - 1; zidx3 = i; zidx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			zidy1 = ny - 1; zidy2 = 0; zidy3 = ny - 1; zidy4 = 0;
		}
		else {
			zidy1 = j - 1; zidy2 = j; zidy3 = j - 1; zidy4 = j;
		}

		zidz1 = k; zidz2 = k; zidz3 = k; zidz4 = k;

		dHydx = (DHy_em_store(zidx3, zidy2, zidz1) - DHy_em_store(zidx1, zidy2, zidz1)) / dx;
		dHxdy = (DHx_em_store(zidx3, zidy2, zidz1) - DHx_em_store(zidx3, zidy1, zidz1)) / dy;

		Jfz = 0.; Jpz = 0.;
		dPz = 0.;
		
		
		
		xcond = 0.;
		ycond = 0.;
		zcond = 0.;

		xer11 = 0.; xer12 = 0.; xer13 = 0.;
		xer21 = 0.; xer22 = 0.; xer23 = 0.;
		xer31 = 0.; xer32 = 0.; xer33 = 0.;

		yer11 = 0.; yer12 = 0.; yer13 = 0.;
		yer21 = 0.; yer22 = 0.; yer23 = 0.;
		yer31 = 0.; yer32 = 0.; yer33 = 0.;

		zer11 = 0.; zer12 = 0.; zer13 = 0.;
		zer21 = 0.; zer22 = 0.; zer23 = 0.;
		zer31 = 0.; zer32 = 0.; zer33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		Px_count = 0.; Py_count = 0.; Pz_count = 0.;



		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(xidx1, xidy1, xidz1);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(xidx1, xidy1, xidz1); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx1, xidy1, xidz1);
		Jpx = Jpx + Jpx_n1_store(xidx1, xidy1, xidz1) + \
			Jpx_n2_store(xidx1, xidy1, xidz1) + Jpx_n3_store(xidx1, xidy1, xidz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx1, xidy1, xidz1);


		mat_type = pt_glb->material_cell(yidx1, yidy1, yidz1);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk2(yidx1, yidy1, yidz1); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx1, yidy1, yidz1);
		Jpy = Jpy + Jpy_n1_store(yidx1, yidy1, yidz1) + \
			Jpy_n2_store(yidx1, yidy1, yidz1) + Jpy_n3_store(yidx1, yidy1, yidz1);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx1, yidy1, yidz1);


		mat_type = pt_glb->material_cell(zidx1, zidy1, zidz1);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk2(zidx1, zidy1, zidz1); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx1, zidy1, zidz1);
		Jpz = Jpz + Jpz_n1_store(zidx1, zidy1, zidz1) + \
			Jpz_n2_store(zidx1, zidy1, zidz1) + Jpz_n3_store(zidx1, zidy1, zidz1);




		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(xidx2, xidy2, xidz2);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(xidx2, xidy2, xidz2); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx2, xidy2, xidz2);
		Jpx = Jpx + Jpx_n1_store(xidx2, xidy2, xidz2) + \
			Jpx_n2_store(xidx2, xidy2, xidz2) + Jpx_n3_store(xidx2, xidy2, xidz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx2, xidy2, xidz2);


		mat_type = pt_glb->material_cell(yidx2, yidy2, yidz2);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk2(yidx2, yidy2, yidz2); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx2, yidy2, yidz2);
		Jpy = Jpy + Jpy_n1_store(yidx2, yidy2, yidz2) + \
			Jpy_n2_store(yidx2, yidy2, yidz2) + Jpy_n3_store(yidx2, yidy2, yidz2);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx2, yidy2, yidz2);


		mat_type = pt_glb->material_cell(zidx2, zidy2, zidz2);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk2(zidx2, zidy2, zidz2); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx2, zidy2, zidz2);
		Jpz = Jpz + Jpz_n1_store(zidx2, zidy2, zidz2) + \
			Jpz_n2_store(zidx2, zidy2, zidz2) + Jpz_n3_store(zidx2, zidy2, zidz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(xidx3, xidy3, xidz3);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(xidx3, xidy3, xidz3); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx3, xidy3, xidz3);
		Jpx = Jpx + Jpx_n1_store(xidx3, xidy3, xidz3) + \
			Jpx_n2_store(xidx3, xidy3, xidz3) + Jpx_n3_store(xidx3, xidy3, xidz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx3, xidy3, xidz3);


		mat_type = pt_glb->material_cell(yidx3, yidy3, yidz3);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk2(yidx3, yidy3, yidz3); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx3, yidy3, yidz3);
		Jpy = Jpy + Jpy_n1_store(yidx3, yidy3, yidz3) + \
			Jpy_n2_store(yidx3, yidy3, yidz3) + Jpy_n3_store(yidx3, yidy3, yidz3);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx3, yidy3, yidz3);


		mat_type = pt_glb->material_cell(zidx3, zidy3, zidz3);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk2(zidx3, zidy3, zidz3); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx3, zidy3, zidz3);
		Jpz = Jpz + Jpz_n1_store(zidx3, zidy3, zidz3) + \
			Jpz_n2_store(zidx3, zidy3, zidz3) + Jpz_n3_store(zidx3, zidy3, zidz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(xidx4, xidy4, xidz4);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk2(xidx4, xidy4, xidz4); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx4, xidy4, xidz4);
		Jpx = Jpx + Jpx_n1_store(xidx4, xidy4, xidz4) + \
			Jpx_n2_store(xidx4, xidy4, xidz4) + Jpx_n3_store(xidx4, xidy4, xidz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx4, xidy4, xidz4);


		mat_type = pt_glb->material_cell(yidx4, yidy4, yidz4);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk2(yidx4, yidy4, yidz4); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx4, yidy4, yidz4);
		Jpy = Jpy + Jpy_n1_store(yidx4, yidy4, yidz4) + \
			Jpy_n2_store(yidx4, yidy4, yidz4) + Jpy_n3_store(yidx4, yidy4, yidz4);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx4, yidy4, yidz4);


		mat_type = pt_glb->material_cell(zidx4, zidy4, zidz4);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk2(zidx4, zidy4, zidz4); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx4, zidy4, zidz4);
		Jpz = Jpz + Jpz_n1_store(zidx4, zidy4, zidz4) + \
			Jpz_n2_store(zidx4, zidy4, zidz4) + Jpz_n3_store(zidx4, zidy4, zidz4);


		if (Px_count < 0.5) {
			Px_count = 1.;
		}
		if (Py_count < 0.5) {
			Py_count = 1.;
		}
		if (Pz_count < 0.5) {
			Pz_count = 1.;
		}

		Jfx = Jfx / 4.;	Jfy = Jfy / 4.;	Jfz = Jfz / 4.;
		Jpx = Jpx / 4.;	Jpy = Jpy / 4.;	Jpz = Jpz / 4.;	
		Jishex = Jishex / 4.;	Jishey = Jishey / 4.;	
		dPx = dPx / Px_count;	dPy = dPy / Py_count;	dPz = dPz / Pz_count;
		
		xcond = xcond / 4.;	ycond = ycond / 4.;	zcond = zcond / 4.;

		xer11 = xer11 / 4.; xer12 = xer12 / 4.; xer13 = xer13 / 4.;
		xer21 = xer21 / 4.; xer22 = xer22 / 4.; xer23 = xer23 / 4.;
		xer31 = xer31 / 4.; xer32 = xer32 / 4.; xer33 = xer33 / 4.;

		yer11 = yer11 / 4.; yer12 = yer12 / 4.; yer13 = yer13 / 4.;
		yer21 = yer21 / 4.; yer22 = yer22 / 4.; yer23 = yer23 / 4.;
		yer31 = yer31 / 4.; yer32 = yer32 / 4.; yer33 = yer33 / 4.;

		zer11 = zer11 / 4.; zer12 = zer12 / 4.; zer13 = zer13 / 4.;
		zer21 = zer21 / 4.; zer22 = zer22 / 4.; zer23 = zer23 / 4.;
		zer31 = zer31 / 4.; zer32 = zer32 / 4.; zer33 = zer33 / 4.;

		xDenominator = -xer13*xer22*xer31 + xer12*xer23*xer31 + xer13*xer21*xer32 -\
			xer11*xer23*xer32 - xer12*xer21*xer33 + xer11*xer22*xer33;
		yDenominator = -yer13*yer22*yer31 + yer12*yer23*yer31 + yer13*yer21*yer32 -\
			yer11*yer23*yer32 - yer12*yer21*yer33 + yer11*yer22*yer33;
		zDenominator = -zer13*zer22*zer31 + zer12*zer23*zer31 + zer13*zer21*zer32 -\
			zer11*zer23*zer32 - zer12*zer21*zer33 + zer11*zer22*zer33;

		xinverse_er11 = (-xer23*xer32 + xer22*xer33) / xDenominator;
		xinverse_er12 = ( xer13*xer32 - xer12*xer33) / xDenominator;
		xinverse_er13 = (-xer13*xer22 + xer12*xer23) / xDenominator;
		xinverse_er21 = ( xer23*xer31 - xer21*xer33) / xDenominator;
		xinverse_er22 = (-xer13*xer31 + xer11*xer33) / xDenominator;
		xinverse_er23 = ( xer13*xer21 - xer11*xer23) / xDenominator;
		xinverse_er31 = (-xer22*xer31 + xer21*xer32) / xDenominator;
		xinverse_er32 = ( xer12*xer31 - xer11*xer32) / xDenominator;
		xinverse_er33 = (-xer12*xer21 + xer11*xer22) / xDenominator;


		yinverse_er11 = (-yer23*yer32 + yer22*yer33) / yDenominator;
		yinverse_er12 = ( yer13*yer32 - yer12*yer33) / yDenominator;
		yinverse_er13 = (-yer13*yer22 + yer12*yer23) / yDenominator;
		yinverse_er21 = ( yer23*yer31 - yer21*yer33) / yDenominator;
		yinverse_er22 = (-yer13*yer31 + yer11*yer33) / yDenominator;
		yinverse_er23 = ( yer13*yer21 - yer11*yer23) / yDenominator;
		yinverse_er31 = (-yer22*yer31 + yer21*yer32) / yDenominator;
		yinverse_er32 = ( yer12*yer31 - yer11*yer32) / yDenominator;
		yinverse_er33 = (-yer12*yer21 + yer11*yer22) / yDenominator;


		zinverse_er11 = (-zer23*zer32 + zer22*zer33) / zDenominator;
		zinverse_er12 = ( zer13*zer32 - zer12*zer33) / zDenominator;
		zinverse_er13 = (-zer13*zer22 + zer12*zer23) / zDenominator;
		zinverse_er21 = ( zer23*zer31 - zer21*zer33) / zDenominator;
		zinverse_er22 = (-zer13*zer31 + zer11*zer33) / zDenominator;
		zinverse_er23 = ( zer13*zer21 - zer11*zer23) / zDenominator;
		zinverse_er31 = (-zer22*zer31 + zer21*zer32) / zDenominator;
		zinverse_er32 = ( zer12*zer31 - zer11*zer32) / zDenominator;
		zinverse_er33 = (-zer12*zer21 + zer11*zer22) / zDenominator;



		dDEx_em_rk2(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er11 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er12 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er13 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEy_em_rk2(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er21 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er22 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er23 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEz_em_rk2(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er31 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er32 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er33 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK
	}
	//-----------END X Y Z component-----------//

}

void EMdynamic_system::get_dE_RK3() {
	// double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	// double dP;

	// long int i, j, k;
	// unsigned int mat_type;
	// material* mat;
	// long int idx1, idy1, idz1;
	// long int idx2, idy2, idz2;
	// long int idx3, idy3, idz3;
	// long int idx4, idy4, idz4;

	// double er;
	// double cond;
	// double Jf, Jp, Jishe;
	// //double Jf_count, Jp_count, Ji_count;
	// double P_count;

	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int xidx1, xidy1, xidz1, yidx1, yidy1, yidz1, zidx1, zidy1, zidz1;
	long int xidx2, xidy2, xidz2, yidx2, yidy2, yidz2, zidx2, zidy2, zidz2;
	long int xidx3, xidy3, xidz3, yidx3, yidy3, yidz3, zidx3, zidy3, zidz3;
	long int xidx4, xidy4, xidz4, yidx4, yidy4, yidz4, zidx4, zidy4, zidz4;

	double xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33;
	double yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33;
	double zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33;

	double xDenominator, yDenominator, zDenominator;

	double xinverse_er11, xinverse_er12, xinverse_er13;
	double xinverse_er21, xinverse_er22, xinverse_er23;
	double xinverse_er31, xinverse_er32, xinverse_er33;
	double yinverse_er11, yinverse_er12, yinverse_er13;  
	double yinverse_er21, yinverse_er22, yinverse_er23;
	double yinverse_er31, yinverse_er32, yinverse_er33;
	double zinverse_er11, zinverse_er12, zinverse_er13;
	double zinverse_er21, zinverse_er22, zinverse_er23;
	double zinverse_er31, zinverse_er32, zinverse_er33;



	double xcond, ycond, zcond;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double Px_count, Py_count, Pz_count;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_half(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK3();
	}

#pragma acc wait
	//-----------X component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k,P_count) default(present) async(1)
// 	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
// 		i = id / ((ny + 1) * (nz + 1));
// 		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

// 		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
// 		}
// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
// 		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk3(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx1, idy1, idz1);
// 		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk3(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx2, idy2, idz2);
// 		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk3(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx3, idy3, idz3);
// 		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk3(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx4, idy4, idz4);
// 		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEx_em_rk3(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END X component-----------//

// 	//-----------Y component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
// 	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
// 		i = id / ((ny) * (nz + 1));
// 		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
// 		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk3(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx1, idy1, idz1);
// 		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk3(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx2, idy2, idz2);
// 		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk3(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx3, idy3, idz3);
// 		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk3(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx4, idy4, idz4);
// 		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEy_em_rk3(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}

// 	//-----------END Y component-----------//

// 	//-----------Z component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
// 	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
// 		i = id / ((ny + 1) * (nz));
// 		j = (id - i * ((ny + 1) * (nz))) / (nz);
// 		k = id - i * ((ny + 1) * (nz)) - j * (nz);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
// 		}

// 		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

// 		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
// 		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

// 		Jf = 0.; Jp = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk3(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx1, idy1, idz1);
// 		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk3(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx2, idy2, idz2);
// 		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk3(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx3, idy3, idz3);
// 		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk3(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx4, idy4, idz4);
// 		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEz_em_rk3(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END Z component-----------//


	//-----------X Y Z component-----------//
#pragma acc parallel loop gang vector \
private(xidx1, xidy1, xidz1, xidx2, xidy2, xidz2, xidx3, xidy3, xidz3, xidx4, xidy4, xidz4,\
yidx1, yidy1, yidz1, yidx2, yidy2, yidz2, yidx3, yidy3, yidz3, yidx4, yidy4, yidz4,\
zidx1, zidy1, zidz1, zidx2, zidy2, zidz2, zidx3, zidy3, zidz3, zidx4, zidy4, zidz4,\
xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33,\
yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33,\
zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33,\
xDenominator,\
xinverse_er11, xinverse_er12, xinverse_er13,\
xinverse_er21, xinverse_er22, xinverse_er23,\
xinverse_er31, xinverse_er32, xinverse_er33,\
yDenominator,\
yinverse_er11, yinverse_er12, yinverse_er13,\
yinverse_er21, yinverse_er22, yinverse_er23,\
yinverse_er31, yinverse_er32, yinverse_er33,\
zDenominator,\
zinverse_er11, zinverse_er12, zinverse_er13,\
zinverse_er21, zinverse_er22, zinverse_er23,\
zinverse_er31, zinverse_er32, zinverse_er33,\
xcond, ycond, zcond,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,\
Px_count,Py_count,Pz_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);
		
		
		// Take the x component
		xidx1 = i; xidx2 = i; xidx3 = i; xidx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			xidy1 = ny - 1; xidy2 = ny - 1; xidy3 = 0; xidy4 = 0;
		}
		else {
			xidy1 = j - 1; xidy2 = j - 1; xidy3 = j; xidy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			xidz1 = nz - 1; xidz2 = 0; xidz3 = nz - 1; xidz4 = 0;
		}
		else {
			xidz1 = k - 1; xidz2 = k; xidz3 = k - 1; xidz4 = k;
		}

		dHzdy = (DHz_em_store(xidx1, xidy3, xidz4) - DHz_em_store(xidx1, xidy1, xidz4)) / dy;
		dHydz = (DHy_em_store(xidx1, xidy3, xidz2) - DHy_em_store(xidx1, xidy3, xidz1)) / dz;

		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.; 

		
		// Take the y component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			yidx1 = nx - 1; yidx2 = nx - 1; yidx3 = 0; yidx4 = 0;
		}
		else {
			yidx1 = i - 1; yidx2 = i - 1; yidx3 = i; yidx4 = i;
		}

		yidy1 = j; yidy2 = j; yidy3 = j; yidy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			yidz1 = nz - 1; yidz2 = 0; yidz3 = nz - 1; yidz4 = 0;
		}
		else {
			yidz1 = k - 1; yidz2 = k; yidz3 = k - 1; yidz4 = k;
		}

		dHxdz = (DHx_em_store(yidx3, yidy1, yidz4) - DHx_em_store(yidx3, yidy1, yidz3)) / dz;
		dHzdx = (DHz_em_store(yidx3, yidy1, yidz2) - DHz_em_store(yidx1, yidy1, yidz2)) / dx;

		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.; 

		
		// Take the z component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			zidx1 = nx - 1; zidx2 = nx - 1; zidx3 = 0; zidx4 = 0;
		}
		else {
			zidx1 = i - 1; zidx2 = i - 1; zidx3 = i; zidx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			zidy1 = ny - 1; zidy2 = 0; zidy3 = ny - 1; zidy4 = 0;
		}
		else {
			zidy1 = j - 1; zidy2 = j; zidy3 = j - 1; zidy4 = j;
		}

		zidz1 = k; zidz2 = k; zidz3 = k; zidz4 = k;

		dHydx = (DHy_em_store(zidx3, zidy2, zidz1) - DHy_em_store(zidx1, zidy2, zidz1)) / dx;
		dHxdy = (DHx_em_store(zidx3, zidy2, zidz1) - DHx_em_store(zidx3, zidy1, zidz1)) / dy;

		Jfz = 0.; Jpz = 0.;
		dPz = 0.;
		
		
		
		xcond = 0.;
		ycond = 0.;
		zcond = 0.;

		xer11 = 0.; xer12 = 0.; xer13 = 0.;
		xer21 = 0.; xer22 = 0.; xer23 = 0.;
		xer31 = 0.; xer32 = 0.; xer33 = 0.;

		yer11 = 0.; yer12 = 0.; yer13 = 0.;
		yer21 = 0.; yer22 = 0.; yer23 = 0.;
		yer31 = 0.; yer32 = 0.; yer33 = 0.;

		zer11 = 0.; zer12 = 0.; zer13 = 0.;
		zer21 = 0.; zer22 = 0.; zer23 = 0.;
		zer31 = 0.; zer32 = 0.; zer33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		Px_count = 0.; Py_count = 0.; Pz_count = 0.;



		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(xidx1, xidy1, xidz1);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(xidx1, xidy1, xidz1); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx1, xidy1, xidz1);
		Jpx = Jpx + Jpx_n1_store(xidx1, xidy1, xidz1) + \
			Jpx_n2_store(xidx1, xidy1, xidz1) + Jpx_n3_store(xidx1, xidy1, xidz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx1, xidy1, xidz1);


		mat_type = pt_glb->material_cell(yidx1, yidy1, yidz1);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk3(yidx1, yidy1, yidz1); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx1, yidy1, yidz1);
		Jpy = Jpy + Jpy_n1_store(yidx1, yidy1, yidz1) + \
			Jpy_n2_store(yidx1, yidy1, yidz1) + Jpy_n3_store(yidx1, yidy1, yidz1);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx1, yidy1, yidz1);


		mat_type = pt_glb->material_cell(zidx1, zidy1, zidz1);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk3(zidx1, zidy1, zidz1); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx1, zidy1, zidz1);
		Jpz = Jpz + Jpz_n1_store(zidx1, zidy1, zidz1) + \
			Jpz_n2_store(zidx1, zidy1, zidz1) + Jpz_n3_store(zidx1, zidy1, zidz1);




		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(xidx2, xidy2, xidz2);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(xidx2, xidy2, xidz2); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx2, xidy2, xidz2);
		Jpx = Jpx + Jpx_n1_store(xidx2, xidy2, xidz2) + \
			Jpx_n2_store(xidx2, xidy2, xidz2) + Jpx_n3_store(xidx2, xidy2, xidz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx2, xidy2, xidz2);


		mat_type = pt_glb->material_cell(yidx2, yidy2, yidz2);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk3(yidx2, yidy2, yidz2); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx2, yidy2, yidz2);
		Jpy = Jpy + Jpy_n1_store(yidx2, yidy2, yidz2) + \
			Jpy_n2_store(yidx2, yidy2, yidz2) + Jpy_n3_store(yidx2, yidy2, yidz2);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx2, yidy2, yidz2);


		mat_type = pt_glb->material_cell(zidx2, zidy2, zidz2);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk3(zidx2, zidy2, zidz2); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx2, zidy2, zidz2);
		Jpz = Jpz + Jpz_n1_store(zidx2, zidy2, zidz2) + \
			Jpz_n2_store(zidx2, zidy2, zidz2) + Jpz_n3_store(zidx2, zidy2, zidz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(xidx3, xidy3, xidz3);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(xidx3, xidy3, xidz3); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx3, xidy3, xidz3);
		Jpx = Jpx + Jpx_n1_store(xidx3, xidy3, xidz3) + \
			Jpx_n2_store(xidx3, xidy3, xidz3) + Jpx_n3_store(xidx3, xidy3, xidz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx3, xidy3, xidz3);


		mat_type = pt_glb->material_cell(yidx3, yidy3, yidz3);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk3(yidx3, yidy3, yidz3); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx3, yidy3, yidz3);
		Jpy = Jpy + Jpy_n1_store(yidx3, yidy3, yidz3) + \
			Jpy_n2_store(yidx3, yidy3, yidz3) + Jpy_n3_store(yidx3, yidy3, yidz3);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx3, yidy3, yidz3);


		mat_type = pt_glb->material_cell(zidx3, zidy3, zidz3);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk3(zidx3, zidy3, zidz3); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx3, zidy3, zidz3);
		Jpz = Jpz + Jpz_n1_store(zidx3, zidy3, zidz3) + \
			Jpz_n2_store(zidx3, zidy3, zidz3) + Jpz_n3_store(zidx3, zidy3, zidz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(xidx4, xidy4, xidz4);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk3(xidx4, xidy4, xidz4); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx4, xidy4, xidz4);
		Jpx = Jpx + Jpx_n1_store(xidx4, xidy4, xidz4) + \
			Jpx_n2_store(xidx4, xidy4, xidz4) + Jpx_n3_store(xidx4, xidy4, xidz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx4, xidy4, xidz4);


		mat_type = pt_glb->material_cell(yidx4, yidy4, yidz4);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk3(yidx4, yidy4, yidz4); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx4, yidy4, yidz4);
		Jpy = Jpy + Jpy_n1_store(yidx4, yidy4, yidz4) + \
			Jpy_n2_store(yidx4, yidy4, yidz4) + Jpy_n3_store(yidx4, yidy4, yidz4);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx4, yidy4, yidz4);


		mat_type = pt_glb->material_cell(zidx4, zidy4, zidz4);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk3(zidx4, zidy4, zidz4); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx4, zidy4, zidz4);
		Jpz = Jpz + Jpz_n1_store(zidx4, zidy4, zidz4) + \
			Jpz_n2_store(zidx4, zidy4, zidz4) + Jpz_n3_store(zidx4, zidy4, zidz4);


		if (Px_count < 0.5) {
			Px_count = 1.;
		}
		if (Py_count < 0.5) {
			Py_count = 1.;
		}
		if (Pz_count < 0.5) {
			Pz_count = 1.;
		}

		Jfx = Jfx / 4.;	Jfy = Jfy / 4.;	Jfz = Jfz / 4.;
		Jpx = Jpx / 4.;	Jpy = Jpy / 4.;	Jpz = Jpz / 4.;	
		Jishex = Jishex / 4.;	Jishey = Jishey / 4.;	
		dPx = dPx / Px_count;	dPy = dPy / Py_count;	dPz = dPz / Pz_count;
		
		xcond = xcond / 4.;	ycond = ycond / 4.;	zcond = zcond / 4.;

		xer11 = xer11 / 4.; xer12 = xer12 / 4.; xer13 = xer13 / 4.;
		xer21 = xer21 / 4.; xer22 = xer22 / 4.; xer23 = xer23 / 4.;
		xer31 = xer31 / 4.; xer32 = xer32 / 4.; xer33 = xer33 / 4.;

		yer11 = yer11 / 4.; yer12 = yer12 / 4.; yer13 = yer13 / 4.;
		yer21 = yer21 / 4.; yer22 = yer22 / 4.; yer23 = yer23 / 4.;
		yer31 = yer31 / 4.; yer32 = yer32 / 4.; yer33 = yer33 / 4.;

		zer11 = zer11 / 4.; zer12 = zer12 / 4.; zer13 = zer13 / 4.;
		zer21 = zer21 / 4.; zer22 = zer22 / 4.; zer23 = zer23 / 4.;
		zer31 = zer31 / 4.; zer32 = zer32 / 4.; zer33 = zer33 / 4.;

		xDenominator = -xer13*xer22*xer31 + xer12*xer23*xer31 + xer13*xer21*xer32 - \
			xer11*xer23*xer32 - xer12*xer21*xer33 + xer11*xer22*xer33;
		yDenominator = -yer13*yer22*yer31 + yer12*yer23*yer31 + yer13*yer21*yer32 - \
			yer11*yer23*yer32 - yer12*yer21*yer33 + yer11*yer22*yer33;
		zDenominator = -zer13*zer22*zer31 + zer12*zer23*zer31 + zer13*zer21*zer32 - \
			zer11*zer23*zer32 - zer12*zer21*zer33 + zer11*zer22*zer33;

		xinverse_er11 = (-xer23*xer32 + xer22*xer33) / xDenominator;
		xinverse_er12 = ( xer13*xer32 - xer12*xer33) / xDenominator;
		xinverse_er13 = (-xer13*xer22 + xer12*xer23) / xDenominator;
		xinverse_er21 = ( xer23*xer31 - xer21*xer33) / xDenominator;
		xinverse_er22 = (-xer13*xer31 + xer11*xer33) / xDenominator;
		xinverse_er23 = ( xer13*xer21 - xer11*xer23) / xDenominator;
		xinverse_er31 = (-xer22*xer31 + xer21*xer32) / xDenominator;
		xinverse_er32 = ( xer12*xer31 - xer11*xer32) / xDenominator;
		xinverse_er33 = (-xer12*xer21 + xer11*xer22) / xDenominator;


		yinverse_er11 = (-yer23*yer32 + yer22*yer33) / yDenominator;
		yinverse_er12 = ( yer13*yer32 - yer12*yer33) / yDenominator;
		yinverse_er13 = (-yer13*yer22 + yer12*yer23) / yDenominator;
		yinverse_er21 = ( yer23*yer31 - yer21*yer33) / yDenominator;
		yinverse_er22 = (-yer13*yer31 + yer11*yer33) / yDenominator;
		yinverse_er23 = ( yer13*yer21 - yer11*yer23) / yDenominator;
		yinverse_er31 = (-yer22*yer31 + yer21*yer32) / yDenominator;
		yinverse_er32 = ( yer12*yer31 - yer11*yer32) / yDenominator;
		yinverse_er33 = (-yer12*yer21 + yer11*yer22) / yDenominator;


		zinverse_er11 = (-zer23*zer32 + zer22*zer33) / zDenominator;
		zinverse_er12 = ( zer13*zer32 - zer12*zer33) / zDenominator;
		zinverse_er13 = (-zer13*zer22 + zer12*zer23) / zDenominator;
		zinverse_er21 = ( zer23*zer31 - zer21*zer33) / zDenominator;
		zinverse_er22 = (-zer13*zer31 + zer11*zer33) / zDenominator;
		zinverse_er23 = ( zer13*zer21 - zer11*zer23) / zDenominator;
		zinverse_er31 = (-zer22*zer31 + zer21*zer32) / zDenominator;
		zinverse_er32 = ( zer12*zer31 - zer11*zer32) / zDenominator;
		zinverse_er33 = (-zer12*zer21 + zer11*zer22) / zDenominator;



		dDEx_em_rk3(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er11 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er12 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er13 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEy_em_rk3(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er21 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex -\
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er22 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er23 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEz_em_rk3(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er31 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex -\
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er32 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er33 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK
	}
	//-----------END X Y Z component-----------//




}

void EMdynamic_system::get_dE_RK4() {
	// double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	// double dP;

	// long int i, j, k;
	// unsigned int mat_type;
	// material* mat;
	// long int idx1, idy1, idz1;
	// long int idx2, idy2, idz2;
	// long int idx3, idy3, idz3;
	// long int idx4, idy4, idz4;

	// double er;
	// double cond;
	// double Jf, Jp, Jishe;
	// //double Jf_count, Jp_count, Ji_count;
	// double P_count;

	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dPx, dPy, dPz;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int xidx1, xidy1, xidz1, yidx1, yidy1, yidz1, zidx1, zidy1, zidz1;
	long int xidx2, xidy2, xidz2, yidx2, yidy2, yidz2, zidx2, zidy2, zidz2;
	long int xidx3, xidy3, xidz3, yidx3, yidy3, yidz3, zidx3, zidy3, zidz3;
	long int xidx4, xidy4, xidz4, yidx4, yidy4, yidz4, zidx4, zidy4, zidz4;

	double xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33;
	double yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33;
	double zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33;

	double xDenominator, yDenominator, zDenominator;

	double xinverse_er11, xinverse_er12, xinverse_er13;
	double xinverse_er21, xinverse_er22, xinverse_er23;
	double xinverse_er31, xinverse_er32, xinverse_er33;
	double yinverse_er11, yinverse_er12, yinverse_er13;  
	double yinverse_er21, yinverse_er22, yinverse_er23;
	double yinverse_er31, yinverse_er32, yinverse_er33;
	double zinverse_er11, zinverse_er12, zinverse_er13;
	double zinverse_er21, zinverse_er22, zinverse_er23;
	double zinverse_er31, zinverse_er32, zinverse_er33;

	double xcond, ycond, zcond;
	double Jfx, Jpx, Jishex, Jfy, Jpy, Jishey, Jfz, Jpz;
	//double Jf_count, Jp_count, Ji_count;
	double Px_count, Py_count, Pz_count;




	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_full(); //RK
	}
#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK4();
	}

#pragma acc wait
// 	//-----------X component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k,P_count) default(present) async(1)
// 	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
// 		i = id / ((ny + 1) * (nz + 1));
// 		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

// 		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
// 		}
// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
// 		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk4(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx1, idy1, idz1);
// 		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk4(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx2, idy2, idz2);
// 		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk4(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx3, idy3, idz3);
// 		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpx_glb_rk4(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfx(idx4, idy4, idz4);
// 		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEx_em_rk4(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END X component-----------//

// 	//-----------Y component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
// 	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
// 		i = id / ((ny) * (nz + 1));
// 		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
// 		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

// 		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
// 			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
// 		}
// 		else {
// 			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
// 		}

// 		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
// 		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

// 		Jf = 0.; Jp = 0.; Jishe = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk4(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx1, idy1, idz1);
// 		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk4(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx2, idy2, idz2);
// 		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk4(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx3, idy3, idz3);
// 		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpy_glb_rk4(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfy(idx4, idy4, idz4);
// 		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
// 		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		Jishe = Jishe / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEy_em_rk4(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}

// 	//-----------END Y component-----------//

// 	//-----------Z component-----------//
// #pragma acc parallel loop gang vector \
// private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
// er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
// 	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
// 		i = id / ((ny + 1) * (nz));
// 		j = (id - i * ((ny + 1) * (nz))) / (nz);
// 		k = id - i * ((ny + 1) * (nz)) - j * (nz);

// 		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
// 			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
// 		}
// 		else {
// 			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
// 		}

// 		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
// 			continue; // should update Ex using Liao's absorbing boundary surface
// 		}
// 		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
// 			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
// 		}
// 		else {
// 			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
// 		}

// 		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

// 		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
// 		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

// 		Jf = 0.; Jp = 0.;
// 		dP = 0.; er = 0.; cond = 0.;
// 		//Jf_count = 0.; Jp_count = 0.; 
// 		P_count = 0.;
// 		//----------Polarization 1------//
// 		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk4(idx1, idy1, idz1); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx1, idy1, idz1);
// 		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

// 		//----------Polarization 2------//
// 		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk4(idx2, idy2, idz2); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx2, idy2, idz2);
// 		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

// 		//----------Polarization 3------//
// 		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk4(idx3, idy3, idz3); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx3, idy3, idz3);
// 		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

// 		//----------Polarization 4------//
// 		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
// 		if (mat_type == 0) {
// 			er = er + 1.;
// 		}
// 		else {
// 			mat = &(pt_glb->material_parameters[(mat_type)-1]);
// 			if (pt_glb->if_EM_fromP == true) {
// 				dP = dP + pt_fe->dpz_glb_rk4(idx4, idy4, idz4); //RK
// 			}
// 			if (mat->if_FE == true) {
// 				P_count = P_count + 1.;
// 			}
// 			er = er + mat->r_permittivity;
// 			cond = cond + mat->conductivity;
// 		}
// 		Jf = Jf + DJfz(idx4, idy4, idz4);
// 		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

// 		if (P_count < 0.5) {
// 			P_count = 1.;
// 		}

// 		Jf = Jf / 4.;
// 		Jp = Jp / 4.;
// 		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

// 		dDEz_em_rk4(i, j, k) = pt_glb->dt / e0 / er * \
// 			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
// 	}
// 	//-----------END Z component-----------//



	//-----------X Y Z component-----------//
#pragma acc parallel loop gang vector \
private(xidx1, xidy1, xidz1, xidx2, xidy2, xidz2, xidx3, xidy3, xidz3, xidx4, xidy4, xidz4,\
yidx1, yidy1, yidz1, yidx2, yidy2, yidz2, yidx3, yidy3, yidz3, yidx4, yidy4, yidz4,\
zidx1, zidy1, zidz1, zidx2, zidy2, zidz2, zidx3, zidy3, zidz3, zidx4, zidy4, zidz4,\
xer11, xer12, xer13, xer21, xer22, xer23, xer31, xer32, xer33,\
yer11, yer12, yer13, yer21, yer22, yer23, yer31, yer32, yer33,\
zer11, zer12, zer13, zer21, zer22, zer23, zer31, zer32, zer33,\
xDenominator,\
xinverse_er11, xinverse_er12, xinverse_er13,\
xinverse_er21, xinverse_er22, xinverse_er23,\
xinverse_er31, xinverse_er32, xinverse_er33,\
yDenominator,\
yinverse_er11, yinverse_er12, yinverse_er13,\
yinverse_er21, yinverse_er22, yinverse_er23,\
yinverse_er31, yinverse_er32, yinverse_er33,\
zDenominator,\
zinverse_er11, zinverse_er12, zinverse_er13,\
zinverse_er21, zinverse_er22, zinverse_er23,\
zinverse_er31, zinverse_er32, zinverse_er33,\
xcond, ycond, zcond,\
Jfx,Jpx,Jishex, Jfy,Jpy,Jishey, Jfz,Jpz,\
dPx,dPy,dPz,\
mat_type,mat,\
dHzdy,dHydz,dHxdz,dHzdx,dHxdy,dHydx,\
i,j,k,\
Px_count,Py_count,Pz_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);
		
		
		// Take the x component
		xidx1 = i; xidx2 = i; xidx3 = i; xidx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			xidy1 = ny - 1; xidy2 = ny - 1; xidy3 = 0; xidy4 = 0;
		}
		else {
			xidy1 = j - 1; xidy2 = j - 1; xidy3 = j; xidy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			xidz1 = nz - 1; xidz2 = 0; xidz3 = nz - 1; xidz4 = 0;
		}
		else {
			xidz1 = k - 1; xidz2 = k; xidz3 = k - 1; xidz4 = k;
		}

		dHzdy = (DHz_em_store(xidx1, xidy3, xidz4) - DHz_em_store(xidx1, xidy1, xidz4)) / dy;
		dHydz = (DHy_em_store(xidx1, xidy3, xidz2) - DHy_em_store(xidx1, xidy3, xidz1)) / dz;

		Jfx = 0.; Jpx = 0.; Jishex = 0.;
		dPx = 0.; 

		
		// Take the y component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			yidx1 = nx - 1; yidx2 = nx - 1; yidx3 = 0; yidx4 = 0;
		}
		else {
			yidx1 = i - 1; yidx2 = i - 1; yidx3 = i; yidx4 = i;
		}

		yidy1 = j; yidy2 = j; yidy3 = j; yidy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			yidz1 = nz - 1; yidz2 = 0; yidz3 = nz - 1; yidz4 = 0;
		}
		else {
			yidz1 = k - 1; yidz2 = k; yidz3 = k - 1; yidz4 = k;
		}

		dHxdz = (DHx_em_store(yidx3, yidy1, yidz4) - DHx_em_store(yidx3, yidy1, yidz3)) / dz;
		dHzdx = (DHz_em_store(yidx3, yidy1, yidz2) - DHz_em_store(yidx1, yidy1, yidz2)) / dx;

		Jfy = 0.; Jpy = 0.; Jishey = 0.;
		dPy = 0.; 

		
		// Take the z component
		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			zidx1 = nx - 1; zidx2 = nx - 1; zidx3 = 0; zidx4 = 0;
		}
		else {
			zidx1 = i - 1; zidx2 = i - 1; zidx3 = i; zidx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			zidy1 = ny - 1; zidy2 = 0; zidy3 = ny - 1; zidy4 = 0;
		}
		else {
			zidy1 = j - 1; zidy2 = j; zidy3 = j - 1; zidy4 = j;
		}

		zidz1 = k; zidz2 = k; zidz3 = k; zidz4 = k;

		dHydx = (DHy_em_store(zidx3, zidy2, zidz1) - DHy_em_store(zidx1, zidy2, zidz1)) / dx;
		dHxdy = (DHx_em_store(zidx3, zidy2, zidz1) - DHx_em_store(zidx3, zidy1, zidz1)) / dy;

		Jfz = 0.; Jpz = 0.;
		dPz = 0.;
		
		
		
		xcond = 0.;
		ycond = 0.;
		zcond = 0.;

		xer11 = 0.; xer12 = 0.; xer13 = 0.;
		xer21 = 0.; xer22 = 0.; xer23 = 0.;
		xer31 = 0.; xer32 = 0.; xer33 = 0.;

		yer11 = 0.; yer12 = 0.; yer13 = 0.;
		yer21 = 0.; yer22 = 0.; yer23 = 0.;
		yer31 = 0.; yer32 = 0.; yer33 = 0.;

		zer11 = 0.; zer12 = 0.; zer13 = 0.;
		zer21 = 0.; zer22 = 0.; zer23 = 0.;
		zer31 = 0.; zer32 = 0.; zer33 = 0.;

		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		Px_count = 0.; Py_count = 0.; Pz_count = 0.;



		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(xidx1, xidy1, xidz1);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(xidx1, xidy1, xidz1); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx1, xidy1, xidz1);
		Jpx = Jpx + Jpx_n1_store(xidx1, xidy1, xidz1) + \
			Jpx_n2_store(xidx1, xidy1, xidz1) + Jpx_n3_store(xidx1, xidy1, xidz1);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx1, xidy1, xidz1);


		mat_type = pt_glb->material_cell(yidx1, yidy1, yidz1);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk4(yidx1, yidy1, yidz1); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx1, yidy1, yidz1);
		Jpy = Jpy + Jpy_n1_store(yidx1, yidy1, yidz1) + \
			Jpy_n2_store(yidx1, yidy1, yidz1) + Jpy_n3_store(yidx1, yidy1, yidz1);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx1, yidy1, yidz1);


		mat_type = pt_glb->material_cell(zidx1, zidy1, zidz1);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk4(zidx1, zidy1, zidz1); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx1, zidy1, zidz1);
		Jpz = Jpz + Jpz_n1_store(zidx1, zidy1, zidz1) + \
			Jpz_n2_store(zidx1, zidy1, zidz1) + Jpz_n3_store(zidx1, zidy1, zidz1);




		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(xidx2, xidy2, xidz2);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(xidx2, xidy2, xidz2); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx2, xidy2, xidz2);
		Jpx = Jpx + Jpx_n1_store(xidx2, xidy2, xidz2) + \
			Jpx_n2_store(xidx2, xidy2, xidz2) + Jpx_n3_store(xidx2, xidy2, xidz2);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx2, xidy2, xidz2);


		mat_type = pt_glb->material_cell(yidx2, yidy2, yidz2);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk4(yidx2, yidy2, yidz2); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx2, yidy2, yidz2);
		Jpy = Jpy + Jpy_n1_store(yidx2, yidy2, yidz2) + \
			Jpy_n2_store(yidx2, yidy2, yidz2) + Jpy_n3_store(yidx2, yidy2, yidz2);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx2, yidy2, yidz2);


		mat_type = pt_glb->material_cell(zidx2, zidy2, zidz2);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk4(zidx2, zidy2, zidz2); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx2, zidy2, zidz2);
		Jpz = Jpz + Jpz_n1_store(zidx2, zidy2, zidz2) + \
			Jpz_n2_store(zidx2, zidy2, zidz2) + Jpz_n3_store(zidx2, zidy2, zidz2);



		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(xidx3, xidy3, xidz3);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(xidx3, xidy3, xidz3); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx3, xidy3, xidz3);
		Jpx = Jpx + Jpx_n1_store(xidx3, xidy3, xidz3) + \
			Jpx_n2_store(xidx3, xidy3, xidz3) + Jpx_n3_store(xidx3, xidy3, xidz3);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx3, xidy3, xidz3);


		mat_type = pt_glb->material_cell(yidx3, yidy3, yidz3);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk4(yidx3, yidy3, yidz3); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx3, yidy3, yidz3);
		Jpy = Jpy + Jpy_n1_store(yidx3, yidy3, yidz3) + \
			Jpy_n2_store(yidx3, yidy3, yidz3) + Jpy_n3_store(yidx3, yidy3, yidz3);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx3, yidy3, yidz3);


		mat_type = pt_glb->material_cell(zidx3, zidy3, zidz3);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk4(zidx3, zidy3, zidz3); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx3, zidy3, zidz3);
		Jpz = Jpz + Jpz_n1_store(zidx3, zidy3, zidz3) + \
			Jpz_n2_store(zidx3, zidy3, zidz3) + Jpz_n3_store(zidx3, zidy3, zidz3);



		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(xidx4, xidy4, xidz4);
		if (mat_type == 0) {
			xer11 = xer11 + 1.; 	xer12 = xer12 + 0.; 	xer13 = xer13 + 0.;
			xer21 = xer21 + 0.; 	xer22 = xer22 + 1.; 	xer23 = xer23 + 0.;
			xer31 = xer31 + 0.; 	xer32 = xer32 + 0.; 	xer33 = xer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPx = dPx + pt_fe->dpx_glb_rk4(xidx4, xidy4, xidz4); //RK
			}
			if (mat->if_FE == true) {
				Px_count = Px_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			xer11 = xer11 + mat->r_permittivity11; xer12 = xer12 + mat->r_permittivity12; xer13 = xer13 + mat->r_permittivity13;
			xer21 = xer21 + mat->r_permittivity21; xer22 = xer22 + mat->r_permittivity22; xer23 = xer23 + mat->r_permittivity23;
			xer31 = xer31 + mat->r_permittivity31; xer32 = xer32 + mat->r_permittivity32; xer33 = xer33 + mat->r_permittivity33;
			xcond = xcond + mat->conductivity;
		}
		Jfx = Jfx + DJfx(xidx4, xidy4, xidz4);
		Jpx = Jpx + Jpx_n1_store(xidx4, xidy4, xidz4) + \
			Jpx_n2_store(xidx4, xidy4, xidz4) + Jpx_n3_store(xidx4, xidy4, xidz4);
		Jishex = Jishex + pt_mag->Jx_ISHE(xidx4, xidy4, xidz4);


		mat_type = pt_glb->material_cell(yidx4, yidy4, yidz4);
		if (mat_type == 0) {
			yer11 = yer11 + 1.; 	yer12 = yer12 + 0.; 	yer13 = yer13 + 0.;
			yer21 = yer21 + 0.; 	yer22 = yer22 + 1.; 	yer23 = yer23 + 0.;
			yer31 = yer31 + 0.; 	yer32 = yer32 + 0.; 	yer33 = yer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPy = dPy + pt_fe->dpy_glb_rk4(yidx4, yidy4, yidz4); //RK
			}
			if (mat->if_FE == true) {
				Py_count = Py_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			yer11 = yer11 + mat->r_permittivity11; yer12 = yer12 + mat->r_permittivity12; yer13 = yer13 + mat->r_permittivity13;
			yer21 = yer21 + mat->r_permittivity21; yer22 = yer22 + mat->r_permittivity22; yer23 = yer23 + mat->r_permittivity23;
			yer31 = yer31 + mat->r_permittivity31; yer32 = yer32 + mat->r_permittivity32; yer33 = yer33 + mat->r_permittivity33;
			ycond = ycond + mat->conductivity;
		}
		Jfy = Jfy + DJfy(yidx4, yidy4, yidz4);
		Jpy = Jpy + Jpy_n1_store(yidx4, yidy4, yidz4) + \
			Jpy_n2_store(yidx4, yidy4, yidz4) + Jpy_n3_store(yidx4, yidy4, yidz4);
		Jishey = Jishey + pt_mag->Jy_ISHE(yidx4, yidy4, yidz4);


		mat_type = pt_glb->material_cell(zidx4, zidy4, zidz4);
		if (mat_type == 0) {
			zer11 = zer11 + 1.; 	zer12 = zer12 + 0.; 	zer13 = zer13 + 0.;
			zer21 = zer21 + 0.; 	zer22 = zer22 + 1.; 	zer23 = zer23 + 0.;
			zer31 = zer31 + 0.; 	zer32 = zer32 + 0.; 	zer33 = zer33 + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dPz = dPz + pt_fe->dpz_glb_rk4(zidx4, zidy4, zidz4); //RK
			}
			if (mat->if_FE == true) {
				Pz_count = Pz_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			zer11 = zer11 + mat->r_permittivity11; zer12 = zer12 + mat->r_permittivity12; zer13 = zer13 + mat->r_permittivity13;
			zer21 = zer21 + mat->r_permittivity21; zer22 = zer22 + mat->r_permittivity22; zer23 = zer23 + mat->r_permittivity23;
			zer31 = zer31 + mat->r_permittivity31; zer32 = zer32 + mat->r_permittivity32; zer33 = zer33 + mat->r_permittivity33;
			zcond = zcond + mat->conductivity;
		}
		Jfz = Jfz + DJfz(zidx4, zidy4, zidz4);
		Jpz = Jpz + Jpz_n1_store(zidx4, zidy4, zidz4) + \
			Jpz_n2_store(zidx4, zidy4, zidz4) + Jpz_n3_store(zidx4, zidy4, zidz4);


		if (Px_count < 0.5) {
			Px_count = 1.;
		}
		if (Py_count < 0.5) {
			Py_count = 1.;
		}
		if (Pz_count < 0.5) {
			Pz_count = 1.;
		}

		Jfx = Jfx / 4.;	Jfy = Jfy / 4.;	Jfz = Jfz / 4.;
		Jpx = Jpx / 4.;	Jpy = Jpy / 4.;	Jpz = Jpz / 4.;	
		Jishex = Jishex / 4.;	Jishey = Jishey / 4.;	
		dPx = dPx / Px_count;	dPy = dPy / Py_count;	dPz = dPz / Pz_count;
		
		xcond = xcond / 4.;	ycond = ycond / 4.;	zcond = zcond / 4.;

		xer11 = xer11 / 4.; xer12 = xer12 / 4.; xer13 = xer13 / 4.;
		xer21 = xer21 / 4.; xer22 = xer22 / 4.; xer23 = xer23 / 4.;
		xer31 = xer31 / 4.; xer32 = xer32 / 4.; xer33 = xer33 / 4.;

		yer11 = yer11 / 4.; yer12 = yer12 / 4.; yer13 = yer13 / 4.;
		yer21 = yer21 / 4.; yer22 = yer22 / 4.; yer23 = yer23 / 4.;
		yer31 = yer31 / 4.; yer32 = yer32 / 4.; yer33 = yer33 / 4.;

		zer11 = zer11 / 4.; zer12 = zer12 / 4.; zer13 = zer13 / 4.;
		zer21 = zer21 / 4.; zer22 = zer22 / 4.; zer23 = zer23 / 4.;
		zer31 = zer31 / 4.; zer32 = zer32 / 4.; zer33 = zer33 / 4.;

		xDenominator = -xer13*xer22*xer31 + xer12*xer23*xer31 + xer13*xer21*xer32 -\
			xer11*xer23*xer32 - xer12*xer21*xer33 + xer11*xer22*xer33;
		yDenominator = -yer13*yer22*yer31 + yer12*yer23*yer31 + yer13*yer21*yer32 -\
			yer11*yer23*yer32 - yer12*yer21*yer33 + yer11*yer22*yer33;
		zDenominator = -zer13*zer22*zer31 + zer12*zer23*zer31 + zer13*zer21*zer32 -\
			zer11*zer23*zer32 - zer12*zer21*zer33 + zer11*zer22*zer33;

		xinverse_er11 = (-xer23*xer32 + xer22*xer33) / xDenominator;
		xinverse_er12 = ( xer13*xer32 - xer12*xer33) / xDenominator;
		xinverse_er13 = (-xer13*xer22 + xer12*xer23) / xDenominator;
		xinverse_er21 = ( xer23*xer31 - xer21*xer33) / xDenominator;
		xinverse_er22 = (-xer13*xer31 + xer11*xer33) / xDenominator;
		xinverse_er23 = ( xer13*xer21 - xer11*xer23) / xDenominator;
		xinverse_er31 = (-xer22*xer31 + xer21*xer32) / xDenominator;
		xinverse_er32 = ( xer12*xer31 - xer11*xer32) / xDenominator;
		xinverse_er33 = (-xer12*xer21 + xer11*xer22) / xDenominator;


		yinverse_er11 = (-yer23*yer32 + yer22*yer33) / yDenominator;
		yinverse_er12 = ( yer13*yer32 - yer12*yer33) / yDenominator;
		yinverse_er13 = (-yer13*yer22 + yer12*yer23) / yDenominator;
		yinverse_er21 = ( yer23*yer31 - yer21*yer33) / yDenominator;
		yinverse_er22 = (-yer13*yer31 + yer11*yer33) / yDenominator;
		yinverse_er23 = ( yer13*yer21 - yer11*yer23) / yDenominator;
		yinverse_er31 = (-yer22*yer31 + yer21*yer32) / yDenominator;
		yinverse_er32 = ( yer12*yer31 - yer11*yer32) / yDenominator;
		yinverse_er33 = (-yer12*yer21 + yer11*yer22) / yDenominator;


		zinverse_er11 = (-zer23*zer32 + zer22*zer33) / zDenominator;
		zinverse_er12 = ( zer13*zer32 - zer12*zer33) / zDenominator;
		zinverse_er13 = (-zer13*zer22 + zer12*zer23) / zDenominator;
		zinverse_er21 = ( zer23*zer31 - zer21*zer33) / zDenominator;
		zinverse_er22 = (-zer13*zer31 + zer11*zer33) / zDenominator;
		zinverse_er23 = ( zer13*zer21 - zer11*zer23) / zDenominator;
		zinverse_er31 = (-zer22*zer31 + zer21*zer32) / zDenominator;
		zinverse_er32 = ( zer12*zer31 - zer11*zer32) / zDenominator;
		zinverse_er33 = (-zer12*zer21 + zer11*zer22) / zDenominator;



		dDEx_em_rk4(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er11 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er12 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er13 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEy_em_rk4(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er21 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er22 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er23 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK

		dDEz_em_rk4(i, j, k) = \
			pt_glb->dt / e0 * xinverse_er31 * \
			(dHzdy - dHydz - Jfx - Jpx - Jishex - \
			xcond * DEx_em_store(i, j, k) - dPx / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er32 * \
			(dHxdz - dHzdx - Jfy - Jpy - Jishey - \
			ycond * DEy_em_store(i, j, k) - dPy / pt_glb->dt) + \
			pt_glb->dt / e0 * xinverse_er33 * \
			(dHydx - dHxdy - Jfz - Jpz - \
			zcond * DEz_em_store(i, j, k) - dPz / pt_glb->dt); //RK
	}
	//-----------END X Y Z component-----------//





}

//void EMdynamic_system::get_dH() {
//	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
//	double dMx1, dMy1, dMz1;
//	double dMx2, dMy2, dMz2;
//	double Ms1, Ms2;
//
//	long int i, j, k;
//	unsigned int mat_type;
//	material* mat;
//	long int idx1, idy1, idz1;
//	long int idx2, idy2, idz2;
//
//	//------------------X component-------------------//
//#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
//Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k) default(present) async(1)
//	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
//		i = id / ((ny) * (nz));
//		j = (id - i * ((ny) * (nz))) / (nz);
//		k = id - i * ((ny) * (nz)) - j * (nz);
//
//		if (i == 0) {
//			if (pt_geo->periodicX == true) {
//				idx1 = nx - 1;
//			}
//			else {
//				idx1 = 0;
//			}
//		}
//		else {
//			idx1 = i - 1;
//		}
//		if (i == nx) {
//			if (pt_geo->periodicX == true) {
//				idx2 = 0;
//			}
//			else {
//				idx2 = nx - 1;
//			}
//		}
//		else {
//			idx2 = i;
//		}
//
//		idy1 = j; idy2 = j;
//		idz1 = k; idz2 = k;
//		dEzdy = (DEz_em(i, j + 1, k) - DEz_em(i, j, k)) / dy;
//		dEydz = (DEy_em(i, j, k + 1) - DEy_em(i, j, k)) / dz;
//
//		//----------Magnetization 1------//
//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//
//		if (mat_type == 0) {
//			Ms1 = 0.;
//			dMx1 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms1 = mat->Ms;
//			dMx1 = pt_mag->dmx_glb(idx1, idy1, idz1) * Ms1;
//		}
//		//----------Magnetization 2------//
//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//
//		if (mat_type == 0) {
//			Ms2 = 0.;
//			dMx2 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms2 = mat->Ms;
//			dMx2 = pt_mag->dmx_glb(idx2, idy2, idz2) * Ms2;
//		}
//		//---------END calculating temporal change of M--------//
//
//		dDHx_em(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) * 0.5;
//		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
//		//	dDHx_em(i, j, k) = 0.;
//		//}
//	}
//
//	//------------------END X component-------------------//
//
//	//------------------Y component-------------------//
//#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
//Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k) default(present) async(2)
//	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
//		i = id / ((ny + 1) * (nz));
//		j = (id - i * ((ny + 1) * (nz))) / (nz);
//		k = id - i * ((ny + 1) * (nz)) - j * (nz);
//
//		idx1 = i; idx2 = i;
//		if (j == 0) {
//			if (pt_geo->periodicY == true) {
//				idy1 = ny - 1;
//			}
//			else {
//				idy1 = 0;
//			}
//		}
//		else {
//			idy1 = j - 1;
//		}
//		if (j == ny) {
//			if (pt_geo->periodicY == true) {
//				idy2 = 0;
//			}
//			else {
//				idy2 = ny - 1;
//			}
//		}
//		else {
//			idy2 = j;
//		}
//		idz1 = k; idz2 = k;
//
//		dExdz = (DEx_em(i, j, k + 1) - DEx_em(i, j, k)) / dz;
//		dEzdx = (DEz_em(i + 1, j, k) - DEz_em(i, j, k)) / dx;
//
//		//----------Magnetization 1------//
//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//
//		if (mat_type == 0) {
//			dMy1 = 0.;
//			Ms1 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms1 = mat->Ms;
//			dMy1 = pt_mag->dmy_glb(idx1, idy1, idz1) * Ms1;
//		}
//		//----------Magnetization 2------//
//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//
//		if (mat_type == 0) {
//			dMy2 = 0.;
//			Ms2 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms2 = mat->Ms;
//			dMy2 = pt_mag->dmy_glb(idx2, idy2, idz2) * Ms2;
//		}
//		//---------END calculating temporal change of M--------//
//
//		dDHy_em(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) * 0.5;
//		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
//		//	dDHy_em(i, j, k) = 0.;
//		//}
//	}
//	//------------------END Y component-------------------//
//
//	//------------------Z component-------------------//
//#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
//Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k) default(present) async(3)
//	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
//		i = id / ((ny) * (nz + 1));
//		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
//		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);
//
//		dEydx = (DEy_em(i + 1, j, k) - DEy_em(i, j, k)) / dx;
//		dExdy = (DEx_em(i, j + 1, k) - DEx_em(i, j, k)) / dy;
//
//		idx1 = i; idx2 = i;
//		idy1 = j; idy2 = j;
//		//----------Magnetization 1------//
//		if (k == 0) {
//			if (pt_geo->periodicZ == true) {
//				idz1 = nz - 1;
//			}
//			else {
//				idz1 = 0;
//			}
//		}
//		else {
//			idz1 = k - 1;
//		}
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//
//		if (mat_type == 0) {
//			dMz1 = 0.;
//			Ms1 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms1 = mat->Ms;
//			dMz1 = pt_mag->dmz_glb(idx1, idy1, idz1) * Ms1;
//		}
//		//----------Magnetization 2------//
//		if (k == nz) {
//			if (pt_geo->periodicZ == true) {
//				idz2 = 0;
//			}
//			else {
//				idz2 = nz - 1;
//			}
//		}
//		else {
//			idz2 = k;
//		}
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//
//		if (mat_type == 0) {
//			dMz2 = 0.;
//			Ms2 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms2 = mat->Ms;
//			dMz2 = pt_mag->dmz_glb(idx2, idy2, idz2) * Ms2;
//		}
//		//---------END calculating temporal change of M--------//
//
//		dDHz_em(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) * 0.5;
//		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
//		//	dDHz_em(i, j, k) = 0.;
//		//}
//	}
//	//------------------END Z component-------------------//
//}

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

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
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
			if (pt_geo->periodicX == true) {
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
		dDHx_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
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
			if (pt_geo->periodicY == true) {
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
		dDHy_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
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
			if (pt_geo->periodicZ == true) {
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

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
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
			if (pt_geo->periodicX == true) {
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
		dDHx_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
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
			if (pt_geo->periodicY == true) {
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
		dDHy_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
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
			if (pt_geo->periodicZ == true) {
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

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
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
			if (pt_geo->periodicX == true) {
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
		dDHx_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
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
			if (pt_geo->periodicY == true) {
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
		dDHy_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
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
			if (pt_geo->periodicZ == true) {
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

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
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
			if (pt_geo->periodicX == true) {
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
		dDHx_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
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
			if (pt_geo->periodicY == true) {
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
		dDHy_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
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
			if (pt_geo->periodicZ == true) {
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
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk1(id) * 0.5;
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
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk2(id) * 0.5;
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
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk3(id);
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
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em(id) = DHy_em(id) + dDHy_em_rk1(id) / 6. + dDHy_em_rk2(id) / 3. + dDHy_em_rk3(id) / 3. + dDHy_em_rk4(id) / 6.;
			DHy_em_store(id) = DHy_em(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em(id) = DHz_em(id) + dDHz_em_rk1(id) / 6. + dDHz_em_rk2(id) / 3. + dDHz_em_rk3(id) / 3. + dDHz_em_rk4(id) / 6.;
			DHz_em_store(id) = DHz_em(id);
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
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk1(id) * 0.5;
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM_half();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
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
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk2(id) * 0.5;
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM_half();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
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
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk3(id);
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM_full();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
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
	if (pt_glb->if_periodic_allsurface == false) {
		transfer_pointer();
	}

#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em(id) = DEx_em(id) + dDEx_em_rk1(id) / 6. + dDEx_em_rk2(id) / 3. + dDEx_em_rk3(id) / 3. + dDEx_em_rk4(id) / 6.;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em(id) = DEy_em(id) + dDEy_em_rk1(id) / 6. + dDEy_em_rk2(id) / 3. + dDEy_em_rk3(id) / 3. + dDEy_em_rk4(id) / 6.;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em(id) = DEz_em(id) + dDEz_em_rk1(id) / 6. + dDEz_em_rk2(id) / 3. + dDEz_em_rk3(id) / 3. + dDEz_em_rk4(id) / 6.;
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
		update_DE_Boundary();
	}

#pragma acc parallel default(present) async(8)
	{
		DEx_em_store = DEx_em;
		DEy_em_store = DEy_em;
		DEz_em_store = DEz_em;
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
					DJfx(i, j, k) = cos(PI * (pt_glb->Jf_rotate_xy +90.)/ 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * (pt_glb->Jf_rotate_xy +90.)/ 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
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
					DJfx(i, j, k) = cos(PI * (pt_glb->Jf_rotate_xy +90.)/ 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * (pt_glb->Jf_rotate_xy +90.)/ 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
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
					DJfx(i, j, k) = cos(PI * (pt_glb->Jf_rotate_xy +90.)/ 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
					DJfy(i, j, k) = sin(PI * (pt_glb->Jf_rotate_xy +90.)/ 180.) * pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sin(2. * PI * pt_glb->Jf_input_freq * temporal_var);
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
