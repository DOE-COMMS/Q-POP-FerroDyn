#include "elastic_system.h"
#include <iostream>

#pragma acc routine gang// nohost
void elastic_system::get_eigenstrain_static() {
	unsigned int mat_type;
	material* mat;
	double mx, my, mz;
	double mx_AFM1, my_AFM1, mz_AFM1;
	double mx_AFM2, my_AFM2, mz_AFM2;
	double px, py, pz;
	double Q11, Q12, Q44;
	double F11, F12, F44;
	double lamda100, lamda111;
	double lamda100_AFM1, lamda111_AFM1;
	double lamda100_AFM2, lamda111_AFM2;

#pragma acc loop gang vector private(mat_type,mx, my, mz,mx_AFM1, my_AFM1, mz_AFM1,mx_AFM2, my_AFM2, mz_AFM2,px, py, pz,\
Q11, Q12, Q44,F11, F12, F44,lamda100, lamda111,lamda100_AFM1, lamda111_AFM1, lamda100_AFM2, lamda111_AFM2,mat)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[mat_type - 1]);
			Q11 = mat->Q11; Q12 = mat->Q12; Q44 = mat->Q44;
			F11 = mat->F11; F12 = mat->F12; F44 = mat->F44;
			lamda100 = mat->lamda100; lamda111 = mat->lamda111;
			lamda100_AFM1 = mat->lamda100_AFM1; lamda111_AFM1 = mat->lamda111_AFM1;
			lamda100_AFM2 = mat->lamda100_AFM2; lamda111_AFM2 = mat->lamda111_AFM2;
			pt_math->transform_vector_glb2crt(pt_fe->px_glb(id), pt_fe->py_glb(id), pt_fe->pz_glb(id), \
				px, py, pz);
			pt_math->transform_vector_glb2crt(pt_mag->mx_glb(id), pt_mag->my_glb(id), pt_mag->mz_glb(id), \
				mx, my, mz);
			pt_math->transform_vector_glb2crt(pt_mag->mx_AFM1_glb(id), pt_mag->my_AFM1_glb(id), pt_mag->mz_AFM1_glb(id), \
				mx_AFM1, my_AFM1, mz_AFM1);
			pt_math->transform_vector_glb2crt(pt_mag->mx_AFM2_glb(id), pt_mag->my_AFM2_glb(id), pt_mag->mz_AFM2_glb(id), \
				mx_AFM2, my_AFM2, mz_AFM2);

			exx0t0_crt(id) = 1.5 * lamda100 * (pow(mx, 2.) - 1. / 3.) + Q11 * pow(px, 2.) + Q12 * (pow(py, 2.) + pow(pz, 2.))\
				+ 0.5 * 1.5 * lamda100_AFM1 * (pow(mx_AFM1, 2.) - 1. / 3.) + 0.5 * 1.5 * lamda100_AFM2 * (pow(mx_AFM2, 2.) - 1. / 3.);

			eyy0t0_crt(id) = 1.5 * lamda100 * (pow(my, 2.) - 1. / 3.) + Q11 * pow(py, 2.) + Q12 * (pow(px, 2.) + pow(pz, 2.))\
				+ 0.5 * 1.5 * lamda100_AFM1 * (pow(my_AFM1, 2.) - 1. / 3.) + 0.5 * 1.5 * lamda100_AFM2 * (pow(my_AFM2, 2.) - 1. / 3.);

			ezz0t0_crt(id) = 1.5 * lamda100 * (pow(mz, 2.) - 1. / 3.) + Q11 * pow(pz, 2.) + Q12 * (pow(py, 2.) + pow(px, 2.))\
				+ 0.5 * 1.5 * lamda100_AFM1 * (pow(mz_AFM1, 2.) - 1. / 3.) + 0.5 * 1.5 * lamda100_AFM2 * (pow(mz_AFM2, 2.) - 1. / 3.);

			eyz0t0_crt(id) = 1.5 * lamda111 * my * mz + Q44 * py * pz \
				+ 0.5 * 1.5 * lamda111_AFM1 * my_AFM1 * mz_AFM1 + 0.5 * 1.5 * lamda111_AFM2 * my_AFM2 * mz_AFM2;

			exz0t0_crt(id) = 1.5 * lamda111 * mx * mz + Q44 * px * pz \
				+ 0.5 * 1.5 * lamda111_AFM1 * mx_AFM1 * mz_AFM1 + 0.5 * 1.5 * lamda111_AFM2 * mx_AFM2 * mz_AFM2;

			exy0t0_crt(id) = 1.5 * lamda111 * my * mx + Q44 * py * px \
				+ 0.5 * 1.5 * lamda111_AFM1 * my_AFM1 * mx_AFM1 + 0.5 * 1.5 * lamda111_AFM2 * my_AFM2 * mx_AFM2;

			if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
				exx0t0_crt(id) = exx0t0_crt(id) \
					- F11 * pt_fe->dpxdx(id) - F12 * (pt_fe->dpzdz(id) + pt_fe->dpydy(id));
				eyy0t0_crt(id) = eyy0t0_crt(id) \
					- F11 * pt_fe->dpydy(id) - F12 * (pt_fe->dpzdz(id) + pt_fe->dpxdx(id));
				ezz0t0_crt(id) = ezz0t0_crt(id) \
					- F11 * pt_fe->dpzdz(id) - F12 * (pt_fe->dpxdx(id) + pt_fe->dpydy(id));

				eyz0t0_crt(id) = eyz0t0_crt(id) \
					- F44 / 2. * (pt_fe->dpydz(id) + pt_fe->dpzdy(id));
				exz0t0_crt(id) = exz0t0_crt(id) \
					- F44 / 2. * (pt_fe->dpxdz(id) + pt_fe->dpzdx(id));
				exy0t0_crt(id) = exy0t0_crt(id) \
					- F44 / 2. * (pt_fe->dpydx(id) + pt_fe->dpxdy(id));
			}
		}
		else {
			exx0t0_crt(id) = 0.;
			eyy0t0_crt(id) = 0.;
			ezz0t0_crt(id) = 0.;
			eyz0t0_crt(id) = 0.;
			exz0t0_crt(id) = 0.;
			exy0t0_crt(id) = 0.;
		}
	}
}

#pragma acc routine gang
void elastic_system::get_eigenstrain_average_film() {
	unsigned long i, j;
	double nxny;
#pragma acc loop gang vector private(i,j,nxny)
	for (unsigned long id = 0; id < nz; id++) {
		exx0_ave_film[id] = 0.;
		eyy0_ave_film[id] = 0.;
		ezz0_ave_film[id] = 0.;
		eyz0_ave_film[id] = 0.;
		exz0_ave_film[id] = 0.;
		exy0_ave_film[id] = 0.;
		for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
				exx0_ave_film[id] += exx0t0_crt(i, j, id);
				eyy0_ave_film[id] += eyy0t0_crt(i, j, id);
				ezz0_ave_film[id] += ezz0t0_crt(i, j, id);
				eyz0_ave_film[id] += eyz0t0_crt(i, j, id);
				exz0_ave_film[id] += exz0t0_crt(i, j, id);
				exy0_ave_film[id] += exy0t0_crt(i, j, id);
			}
		}
		nxny = static_cast<double>(nx * ny);
		exx0_ave_film[id] = exx0_ave_film[id] / nxny;
		eyy0_ave_film[id] = eyy0_ave_film[id] / nxny;
		ezz0_ave_film[id] = ezz0_ave_film[id] / nxny;
		eyz0_ave_film[id] = eyz0_ave_film[id] / nxny;
		exz0_ave_film[id] = exz0_ave_film[id] / nxny;
		exy0_ave_film[id] = exy0_ave_film[id] / nxny;
	}
}

//#pragma acc routine nohost
//void elastic_system::get_Eelas_sum(/*matrix3d<double>& exx0, matrix3d<double>& eyy0, matrix3d<double>& ezz0, \
//	matrix3d<double>& exy0, matrix3d<double>& exz0, matrix3d<double>& eyz0, \
//	matrix3d<double>& exx, matrix3d<double>& eyy, matrix3d<double>& ezz, \
//	matrix3d<double>& exy, matrix3d<double>& exz, matrix3d<double>& eyz, */double& e_sum) {
//
//	unsigned int mat_type;
//	material* mat;
//	double c11, c12, c44;
//	double sxx, syy, szz, sxy, sxz, syz;
//	double e_sum_local;
//
//#pragma acc declare device_resident(mat_type,c11, c12, c44,sxx, syy, szz, sxy, sxz, syz,e_sum_local) deviceptr(mat)
//
//#pragma acc loop seq 
//	for (long int i = 0; i < 1; i++) {
//		e_sum_local = 0.;
//	}
//#pragma acc wait
//
//#pragma acc loop collapse(3) independent private(mat_type,c11, c12, c44,sxx, syy, szz, sxy, sxz, syz) reduction(+:e_sum_local)
//	for (long int i = 0; i < nx; i++) {
//		for (long int j = 0; j < ny; j++) {
//			for (long int k = 0; k < nz; k++) {
//				mat_type = pt_glb->material_cell(i, j, k);
//				if (mat_type != 0) {
//					mat = &(pt_glb->material_parameters[static_cast<long int>(mat_type) - 1]);
//					c11 = mat->c11; c12 = mat->c12; c44 = mat->c44;
//
//					sxx = exxt0_crt(i, j, k) - exx0t0_crt(i, j, k);
//					syy = eyyt0_crt(i, j, k) - eyy0t0_crt(i, j, k);
//					szz = ezzt0_crt(i, j, k) - ezz0t0_crt(i, j, k);
//					syz = eyzt0_crt(i, j, k) - eyz0t0_crt(i, j, k);
//					sxz = exzt0_crt(i, j, k) - exz0t0_crt(i, j, k);
//					sxy = exyt0_crt(i, j, k) - exy0t0_crt(i, j, k);
//
//					e_sum_local = e_sum_local + \
//						0.5 * c11 * (sxx * sxx + syy * syy + szz * szz) \
//						+ c12 * (sxx * syy + sxx * szz + syy * szz) \
//						+ 2. * c44 * (syz * syz + sxz * sxz + sxy * sxy);
//				}
//			}
//		}
//	}
//#pragma acc wait
//#pragma acc loop seq 
//	for (long int i = 0; i < 1; i++) {
//		e_sum = e_sum_local;
//	}
//#pragma acc wait
//}

void elastic_system::get_strain_static() {
	long i, j, k;
	unsigned int itrcount;
	unsigned int mat_type;
	material* mat;
	double c11, c12, c44;

	double sxx, syy, szz, sxy, sxz, syz;

	std::complex<double> uxk, uyk, uzk;
	std::complex<double> stressk_kq[3];
	double kqx, kqy, kqz;

	double e_sum_local;

#pragma acc parallel default(present)
	{
		get_eigenstrain_static();
	}

	if (pt_glb->if_elastostatic == true) {
		if (pt_glb->if_elastostatic_1D == true) {
#pragma acc parallel default(present)
			{
				get_strain_static_1D();
			}
		}
		else {
//			if (pt_glb->if_elastostatic_film == true) {
//#pragma acc parallel default(present)
//				{
//					get_eigenstrain_average_film();
//				}
//			}
#pragma acc parallel default(present)
			{
#pragma acc loop gang vector private(mat_type, mat,c11,c12,c44,i,j,k)
				for (long int id = 0; id < n; id++) {
					if (pt_glb->if_elastostatic_film == true) {
						i = id / (ny * nz);
						j = (id - i * (ny * nz)) / nz;
						k = id - i * (ny * nz) - j * nz;
						//stress0homo_xx(id) = (c11_homo) * (exx0t0_crt(id) - exx0_ave_film[k]) + (c12_homo) * (eyy0t0_crt(id) - eyy0_ave_film[k] + ezz0t0_crt(id) - ezz0_ave_film[k]);
						//stress0homo_yy(id) = (c11_homo) * (eyy0t0_crt(id) - eyy0_ave_film[k]) + (c12_homo) * (exx0t0_crt(id) - exx0_ave_film[k] + ezz0t0_crt(id) - ezz0_ave_film[k]);
						//stress0homo_zz(id) = (c11_homo) * (ezz0t0_crt(id) - ezz0_ave_film[k]) + (c12_homo) * (exx0t0_crt(id) - exx0_ave_film[k] + eyy0t0_crt(id) - eyy0_ave_film[k]);
						//stress0homo_xy(id) = 2. * (c44_homo) * (exy0t0_crt(id) - exy0_ave_film[k]);
						//stress0homo_xz(id) = 2. * (c44_homo) * (exz0t0_crt(id) - exz0_ave_film[k]);
						//stress0homo_yz(id) = 2. * (c44_homo) * (eyz0t0_crt(id) - eyz0_ave_film[k]);
						stress0homo_xx(id) = c11_homo * exx0t0_crt(id) + c12_homo * (eyy0t0_crt(id) + ezz0t0_crt(id));
						stress0homo_yy(id) = c11_homo * eyy0t0_crt(id) + c12_homo * (exx0t0_crt(id) + ezz0t0_crt(id));
						stress0homo_zz(id) = c11_homo * ezz0t0_crt(id) + c12_homo * (exx0t0_crt(id) + eyy0t0_crt(id));
						stress0homo_xy(id) = 2. * c44_homo * exy0t0_crt(id);
						stress0homo_xz(id) = 2. * c44_homo * exz0t0_crt(id);
						stress0homo_yz(id) = 2. * c44_homo * eyz0t0_crt(id);
					}
					else {
						stress0homo_xx(id) = (c11_homo) * (exx0t0_crt(id) + exx_ext_crt(id)) + (c12_homo) * (eyy0t0_crt(id) + eyy_ext_crt(id) + ezz0t0_crt(id) + ezz_ext_crt(id));
						stress0homo_yy(id) = (c11_homo) * (eyy0t0_crt(id) + eyy_ext_crt(id)) + (c12_homo) * (exx0t0_crt(id) + exx_ext_crt(id) + ezz0t0_crt(id) + ezz_ext_crt(id));
						stress0homo_zz(id) = (c11_homo) * (ezz0t0_crt(id) + ezz_ext_crt(id)) + (c12_homo) * (exx0t0_crt(id) + exx_ext_crt(id) + eyy0t0_crt(id) + eyy_ext_crt(id));
						stress0homo_xy(id) = 2. * (c44_homo) * (exy0t0_crt(id) + exy_ext_crt(id));
						stress0homo_xz(id) = 2. * (c44_homo) * (exz0t0_crt(id) + exz_ext_crt(id));
						stress0homo_yz(id) = 2. * (c44_homo) * (eyz0t0_crt(id) + eyz_ext_crt(id));
					}

					mat_type = pt_glb->material_cell(id);
					if (mat_type != 0) {
						mat = &(pt_glb->material_parameters[mat_type - 1]);
						c11 = mat->c11; c12 = mat->c12; c44 = mat->c44;

						if (pt_glb->if_elastostatic_film == true) {
							//stress0_xx(id) = c11 * (exx0t0_crt(id) - exx0_ave_film[k]) + c12 * (eyy0t0_crt(id) - eyy0_ave_film[k] + ezz0t0_crt(id) - ezz0_ave_film[k]);
							//stress0_yy(id) = c11 * (eyy0t0_crt(id) - eyy0_ave_film[k]) + c12 * (exx0t0_crt(id) - exx0_ave_film[k] + ezz0t0_crt(id) - ezz0_ave_film[k]);
							//stress0_zz(id) = c11 * (ezz0t0_crt(id) - ezz0_ave_film[k]) + c12 * (exx0t0_crt(id) - exx0_ave_film[k] + eyy0t0_crt(id) - eyy0_ave_film[k]);
							//stress0_xy(id) = 2. * c44 * (exy0t0_crt(id) - exy0_ave_film[k]);
							//stress0_xz(id) = 2. * c44 * (exz0t0_crt(id) - exz0_ave_film[k]);
							//stress0_yz(id) = 2. * c44 * (eyz0t0_crt(id) - eyz0_ave_film[k]);
							stress0_xx(id) = c11 * exx0t0_crt(id) + c12 * (eyy0t0_crt(id) + ezz0t0_crt(id));
							stress0_yy(id) = c11 * eyy0t0_crt(id) + c12 * (exx0t0_crt(id) + ezz0t0_crt(id));
							stress0_zz(id) = c11 * ezz0t0_crt(id) + c12 * (exx0t0_crt(id) + eyy0t0_crt(id));
							stress0_xy(id) = 2. * c44 * exy0t0_crt(id);
							stress0_xz(id) = 2. * c44 * exz0t0_crt(id);
							stress0_yz(id) = 2. * c44 * eyz0t0_crt(id);
						}
						else {
							stress0_xx(id) = c11 * (exx0t0_crt(id) + exx_ext_crt(id)) + c12 * (eyy0t0_crt(id) + eyy_ext_crt(id) + ezz0t0_crt(id) + ezz_ext_crt(id));
							stress0_yy(id) = c11 * (eyy0t0_crt(id) + eyy_ext_crt(id)) + c12 * (exx0t0_crt(id) + exx_ext_crt(id) + ezz0t0_crt(id) + ezz_ext_crt(id));
							stress0_zz(id) = c11 * (ezz0t0_crt(id) + ezz_ext_crt(id)) + c12 * (exx0t0_crt(id) + exx_ext_crt(id) + eyy0t0_crt(id) + eyy_ext_crt(id));
							stress0_xy(id) = 2. * c44 * (exy0t0_crt(id) + exy_ext_crt(id));
							stress0_xz(id) = 2. * c44 * (exz0t0_crt(id) + exz_ext_crt(id));
							stress0_yz(id) = 2. * c44 * (eyz0t0_crt(id) + eyz_ext_crt(id));

						}
					}
					else {
						stress0_xx(id) = 0.;
						stress0_yy(id) = 0.;
						stress0_zz(id) = 0.;
						stress0_xy(id) = 0.;
						stress0_xz(id) = 0.;
						stress0_yz(id) = 0.;
					}

				}
			}

#pragma acc host_data use_device(	stress0homo_xx.matrix,	stress0homok_xx.matrix,\
									stress0homo_yy.matrix,	stress0homok_yy.matrix,\
									stress0homo_zz.matrix,	stress0homok_zz.matrix,\
									stress0homo_yz.matrix,	stress0homok_yz.matrix,\
									stress0homo_xz.matrix,	stress0homok_xz.matrix,\
									stress0homo_xy.matrix,	stress0homok_xy.matrix)
			{
				cufftExecD2Z(plan_elasto_D2Z, stress0homo_xx.matrix, (cufftDoubleComplex*)(stress0homok_xx.matrix));
				cufftExecD2Z(plan_elasto_D2Z, stress0homo_yy.matrix, (cufftDoubleComplex*)(stress0homok_yy.matrix));
				cufftExecD2Z(plan_elasto_D2Z, stress0homo_zz.matrix, (cufftDoubleComplex*)(stress0homok_zz.matrix));
				cufftExecD2Z(plan_elasto_D2Z, stress0homo_yz.matrix, (cufftDoubleComplex*)(stress0homok_yz.matrix));
				cufftExecD2Z(plan_elasto_D2Z, stress0homo_xz.matrix, (cufftDoubleComplex*)(stress0homok_xz.matrix));
				cufftExecD2Z(plan_elasto_D2Z, stress0homo_xy.matrix, (cufftDoubleComplex*)(stress0homok_xy.matrix));
				//pt_math->fourier_transform3D(stress0homo_xx, stress0homok_xx);
				//pt_math->fourier_transform3D(stress0homo_yy, stress0homok_yy);
				//pt_math->fourier_transform3D(stress0homo_zz, stress0homok_zz);
				//pt_math->fourier_transform3D(stress0homo_yz, stress0homok_yz);
				//pt_math->fourier_transform3D(stress0homo_xz, stress0homok_xz);
				//pt_math->fourier_transform3D(stress0homo_xy, stress0homok_xy);
			}
#pragma acc wait

#pragma acc parallel default(present)
			{
#pragma acc loop gang vector private(i,j,k,kqx,kqy,kqz,stressk_kq,uxk,uyk,uzk)
				for (long int id = 0; id < nx * ny * nz21; id++) {
					i = id / (ny * nz21);
					j = (id - i * (ny * nz21)) / nz21;
					k = id - i * (ny * nz21) - j * nz21;

					if (pt_glb->k_norm(i, j, k) > 1.e-8) {
						pt_math->transform_vector_glb2crt(pt_glb->kx[i], pt_glb->ky[j], pt_glb->kz[k], \
							kqx, kqy, kqz);

						stressk_kq[0] = stress0homok_xx(i, j, k) * kqx + stress0homok_xy(i, j, k) * kqy + stress0homok_xz(i, j, k) * kqz;
						stressk_kq[1] = stress0homok_xy(i, j, k) * kqx + stress0homok_yy(i, j, k) * kqy + stress0homok_yz(i, j, k) * kqz;
						stressk_kq[2] = stress0homok_xz(i, j, k) * kqx + stress0homok_yz(i, j, k) * kqy + stress0homok_zz(i, j, k) * kqz;

						uxk = -im_unit * (Gik_11(i, j, k) * stressk_kq[0] + Gik_12(i, j, k) * stressk_kq[1] + Gik_13(i, j, k) * stressk_kq[2]);
						uyk = -im_unit * (Gik_21(i, j, k) * stressk_kq[0] + Gik_22(i, j, k) * stressk_kq[1] + Gik_23(i, j, k) * stressk_kq[2]);
						uzk = -im_unit * (Gik_31(i, j, k) * stressk_kq[0] + Gik_32(i, j, k) * stressk_kq[1] + Gik_33(i, j, k) * stressk_kq[2]);
					}
					else {
						uxk = 0.; uyk = 0.; uzk = 0.;
					}

					duxdxk(i, j, k) = im_unit * kqx * uxk;
					duxdyk(i, j, k) = im_unit * kqy * uxk;
					duxdzk(i, j, k) = im_unit * kqz * uxk;

					duydxk(i, j, k) = im_unit * kqx * uyk;
					duydyk(i, j, k) = im_unit * kqy * uyk;
					duydzk(i, j, k) = im_unit * kqz * uyk;

					duzdxk(i, j, k) = im_unit * kqx * uzk;
					duzdyk(i, j, k) = im_unit * kqy * uzk;
					duzdzk(i, j, k) = im_unit * kqz * uzk;
				}
			}

#pragma acc host_data use_device(	duxdxk.matrix,	duxdx.matrix,\
									duxdyk.matrix,	duxdy.matrix,\
									duxdzk.matrix,	duxdz.matrix,\
									duydxk.matrix,	duydx.matrix,\
									duydyk.matrix,	duydy.matrix,\
									duydzk.matrix,	duydz.matrix,\
									duzdxk.matrix,	duzdx.matrix,\
									duzdyk.matrix,	duzdy.matrix,\
									duzdzk.matrix,	duzdz.matrix)
			{
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duxdxk.matrix), duxdx.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duxdyk.matrix), duxdy.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duxdzk.matrix), duxdz.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duydxk.matrix), duydx.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duydyk.matrix), duydy.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duydzk.matrix), duydz.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duzdxk.matrix), duzdx.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duzdyk.matrix), duzdy.matrix);
				cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duzdzk.matrix), duzdz.matrix);
				//pt_math->fourier_transform3D(duxdxk, duxdx); duxdx *= scalenn;
				//pt_math->fourier_transform3D(duxdyk, duxdy); duxdy *= scalenn;
				//pt_math->fourier_transform3D(duxdzk, duxdz); duxdz *= scalenn;
				//pt_math->fourier_transform3D(duydxk, duydx); duydx *= scalenn;
				//pt_math->fourier_transform3D(duydyk, duydy); duydy *= scalenn;
				//pt_math->fourier_transform3D(duydzk, duydz); duydz *= scalenn;
				//pt_math->fourier_transform3D(duzdxk, duzdx); duzdx *= scalenn;
				//pt_math->fourier_transform3D(duzdyk, duzdy); duzdy *= scalenn;
				//pt_math->fourier_transform3D(duzdzk, duzdz); duzdz *= scalenn;
			}
#pragma acc wait

#pragma acc parallel default(present)
			{
				duxdx *= scalenn;
				duxdy *= scalenn;
				duxdz *= scalenn;
				duydx *= scalenn;
				duydy *= scalenn;
				duydz *= scalenn;
				duzdx *= scalenn;
				duzdy *= scalenn;
				duzdz *= scalenn;
			}

#pragma acc parallel default(present)
			{
#pragma acc loop gang vector
				for (long int id = 0; id < n; id++) {

					exxt0_crt(id) = duxdx(id);
					eyyt0_crt(id) = duydy(id);
					ezzt0_crt(id) = duzdz(id);
					eyzt0_crt(id) = (duydz(id) + duzdy(id)) / 2.;
					exzt0_crt(id) = (duxdz(id) + duzdx(id)) / 2.;
					exyt0_crt(id) = (duxdy(id) + duydx(id)) / 2.;

				}
			}

			//#pragma acc serial default(present)
			//		{
			//			e_sum = 0.;
			//		}
			//#pragma acc wait
			e_sum_local = 0.;
#pragma acc parallel default(present)
			{
#pragma acc loop gang vector private(mat_type,c11, c12, c44,sxx, syy, szz, sxy, sxz, syz) reduction(+:e_sum_local)
				for (long int id = 0; id < n; id++) {

					mat_type = pt_glb->material_cell(id);
					if (mat_type != 0) {
						mat = &(pt_glb->material_parameters[mat_type - 1]);
						c11 = mat->c11; c12 = mat->c12; c44 = mat->c44;

						sxx = exxt0_crt(id) - exx0t0_crt(id);
						syy = eyyt0_crt(id) - eyy0t0_crt(id);
						szz = ezzt0_crt(id) - ezz0t0_crt(id);
						syz = eyzt0_crt(id) - eyz0t0_crt(id);
						sxz = exzt0_crt(id) - exz0t0_crt(id);
						sxy = exyt0_crt(id) - exy0t0_crt(id);

						e_sum_local = e_sum_local + \
							0.5 * c11 * (sxx * sxx + syy * syy + szz * szz) \
							+ c12 * (sxx * syy + sxx * szz + syy * szz) \
							+ 2. * c44 * (syz * syz + sxz * sxz + sxy * sxy);
					}

				}
				//	get_Eelas_sum(/*exx0t0_crt, eyy0t0_crt, ezz0t0_crt, exy0t0_crt, exz0t0_crt, eyz0t0_crt, exxt0_crt, eyyt0_crt, ezzt0_crt, exyt0_crt, exzt0_crt, eyzt0_crt,*/e_sum);
			}
			e_sum = e_sum_local;

			//#pragma acc serial default(present)
			//		{
			e_elas[0] = 0.; e_elas[1] = 0.; e_elas[2] = 0.; e_elas[3] = 0.;
			e_elas[3] = e_sum;
			//		}
			//#pragma acc wait

			for (itrcount = 0; itrcount < pt_glb->elasto_solver_limit; itrcount++) {
#pragma acc parallel default(present)
				{
#pragma acc loop gang vector private(mat_type, mat,c11,c12,c44)
					for (long int id = 0; id < n; id++) {

						mat_type = pt_glb->material_cell(id);
						if (mat_type != 0) {
							mat = &(pt_glb->material_parameters[mat_type - 1]);
							c11 = mat->c11 - c11_homo; c12 = mat->c12 - c12_homo; c44 = mat->c44 - c44_homo;
						}
						else {
							c11 = 0. - c11_homo; c12 = 0. - c12_homo; c44 = 0. - c44_homo;
						}
						stress_xx(id) = stress0_xx(id) - c11 * exxt0_crt(id) - c12 * (eyyt0_crt(id) + ezzt0_crt(id));
						stress_yy(id) = stress0_yy(id) - c11 * eyyt0_crt(id) - c12 * (exxt0_crt(id) + ezzt0_crt(id));
						stress_zz(id) = stress0_zz(id) - c11 * ezzt0_crt(id) - c12 * (exxt0_crt(id) + eyyt0_crt(id));
						stress_xy(id) = stress0_xy(id) - 2. * c44 * exyt0_crt(id);
						stress_xz(id) = stress0_xz(id) - 2. * c44 * exzt0_crt(id);
						stress_yz(id) = stress0_yz(id) - 2. * c44 * eyzt0_crt(id);

					}
				}

#pragma acc host_data use_device(	stress_xx.matrix,	stressk_xx.matrix,\
									stress_yy.matrix,	stressk_yy.matrix,\
									stress_zz.matrix,	stressk_zz.matrix,\
									stress_yz.matrix,	stressk_yz.matrix,\
									stress_xz.matrix,	stressk_xz.matrix,\
									stress_xy.matrix,	stressk_xy.matrix)
				{
					cufftExecD2Z(plan_elasto_D2Z, stress_xx.matrix, (cufftDoubleComplex*)(stressk_xx.matrix));
					cufftExecD2Z(plan_elasto_D2Z, stress_yy.matrix, (cufftDoubleComplex*)(stressk_yy.matrix));
					cufftExecD2Z(plan_elasto_D2Z, stress_zz.matrix, (cufftDoubleComplex*)(stressk_zz.matrix));
					cufftExecD2Z(plan_elasto_D2Z, stress_yz.matrix, (cufftDoubleComplex*)(stressk_yz.matrix));
					cufftExecD2Z(plan_elasto_D2Z, stress_xz.matrix, (cufftDoubleComplex*)(stressk_xz.matrix));
					cufftExecD2Z(plan_elasto_D2Z, stress_xy.matrix, (cufftDoubleComplex*)(stressk_xy.matrix));

					//pt_math->fourier_transform3D(stress_xx, stressk_xx);
					//pt_math->fourier_transform3D(stress_yy, stressk_yy);
					//pt_math->fourier_transform3D(stress_zz, stressk_zz);
					//pt_math->fourier_transform3D(stress_yz, stressk_yz);
					//pt_math->fourier_transform3D(stress_xz, stressk_xz);
					//pt_math->fourier_transform3D(stress_xy, stressk_xy);
				}
#pragma acc wait

#pragma acc parallel default(present)
				{
#pragma acc loop gang vector private(i,j,k,kqx,kqy,kqz,stressk_kq,uxk,uyk,uzk)
					for (long int id = 0; id < nx * ny * nz21; id++) {
						i = id / (ny * nz21);
						j = (id - i * (ny * nz21)) / nz21;
						k = id - i * (ny * nz21) - j * nz21;

						if (pt_glb->k_norm(i, j, k) > 1.e-8) {
							pt_math->transform_vector_glb2crt(pt_glb->kx[i], pt_glb->ky[j], pt_glb->kz[k], \
								kqx, kqy, kqz);

							stressk_kq[0] = stressk_xx(i, j, k) * kqx + stressk_xy(i, j, k) * kqy + stressk_xz(i, j, k) * kqz;
							stressk_kq[1] = stressk_xy(i, j, k) * kqx + stressk_yy(i, j, k) * kqy + stressk_yz(i, j, k) * kqz;
							stressk_kq[2] = stressk_xz(i, j, k) * kqx + stressk_yz(i, j, k) * kqy + stressk_zz(i, j, k) * kqz;

							uxk = -im_unit * (Gik_11(i, j, k) * stressk_kq[0] + Gik_12(i, j, k) * stressk_kq[1] + Gik_13(i, j, k) * stressk_kq[2]);
							uyk = -im_unit * (Gik_21(i, j, k) * stressk_kq[0] + Gik_22(i, j, k) * stressk_kq[1] + Gik_23(i, j, k) * stressk_kq[2]);
							uzk = -im_unit * (Gik_31(i, j, k) * stressk_kq[0] + Gik_32(i, j, k) * stressk_kq[1] + Gik_33(i, j, k) * stressk_kq[2]);
						}
						else {
							uxk = 0.; uyk = 0.; uzk = 0.;
						}

						duxdxk(i, j, k) = im_unit * kqx * uxk;
						duxdyk(i, j, k) = im_unit * kqy * uxk;
						duxdzk(i, j, k) = im_unit * kqz * uxk;

						duydxk(i, j, k) = im_unit * kqx * uyk;
						duydyk(i, j, k) = im_unit * kqy * uyk;
						duydzk(i, j, k) = im_unit * kqz * uyk;

						duzdxk(i, j, k) = im_unit * kqx * uzk;
						duzdyk(i, j, k) = im_unit * kqy * uzk;
						duzdzk(i, j, k) = im_unit * kqz * uzk;
					}
				}

#pragma acc host_data use_device(	duxdxk.matrix,	duxdx.matrix,\
									duxdyk.matrix,	duxdy.matrix,\
									duxdzk.matrix,	duxdz.matrix,\
									duydxk.matrix,	duydx.matrix,\
									duydyk.matrix,	duydy.matrix,\
									duydzk.matrix,	duydz.matrix,\
									duzdxk.matrix,	duzdx.matrix,\
									duzdyk.matrix,	duzdy.matrix,\
									duzdzk.matrix,	duzdz.matrix)
				{
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duxdxk.matrix), duxdx.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duxdyk.matrix), duxdy.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duxdzk.matrix), duxdz.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duydxk.matrix), duydx.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duydyk.matrix), duydy.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duydzk.matrix), duydz.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duzdxk.matrix), duzdx.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duzdyk.matrix), duzdy.matrix);
					cufftExecZ2D(plan_elasto_Z2D, (cufftDoubleComplex*)(duzdzk.matrix), duzdz.matrix);
					//		pt_math->fourier_transform3D(duxdxk, duxdx); duxdx *= scalenn;
					//		pt_math->fourier_transform3D(duxdyk, duxdy); duxdy *= scalenn;
					//		pt_math->fourier_transform3D(duxdzk, duxdz); duxdz *= scalenn;
					//		pt_math->fourier_transform3D(duydxk, duydx); duydx *= scalenn;
					//		pt_math->fourier_transform3D(duydyk, duydy); duydy *= scalenn;
					//		pt_math->fourier_transform3D(duydzk, duydz); duydz *= scalenn;
					//		pt_math->fourier_transform3D(duzdxk, duzdx); duzdx *= scalenn;
					//		pt_math->fourier_transform3D(duzdyk, duzdy); duzdy *= scalenn;
					//		pt_math->fourier_transform3D(duzdzk, duzdz); duzdz *= scalenn;
				}
#pragma acc wait

#pragma acc parallel default(present)
				{
					duxdx *= scalenn;
					duxdy *= scalenn;
					duxdz *= scalenn;
					duydx *= scalenn;
					duydy *= scalenn;
					duydz *= scalenn;
					duzdx *= scalenn;
					duzdy *= scalenn;
					duzdz *= scalenn;
				}

#pragma acc parallel default(present)
				{
#pragma acc loop gang vector
					for (long int id = 0; id < n; id++) {
						exxt0_crt(id) = duxdx(id);
						eyyt0_crt(id) = duydy(id);
						ezzt0_crt(id) = duzdz(id);
						eyzt0_crt(id) = (duydz(id) + duzdy(id)) / 2.;
						exzt0_crt(id) = (duxdz(id) + duzdx(id)) / 2.;
						exyt0_crt(id) = (duxdy(id) + duydx(id)) / 2.;
					}
				}

				e_sum_local = 0.;
#pragma acc parallel default(present)
				{
#pragma acc loop gang vector private(mat_type,c11, c12, c44,sxx, syy, szz, sxy, sxz, syz) reduction(+:e_sum_local)
					for (long int id = 0; id < n; id++) {
						mat_type = pt_glb->material_cell(id);
						if (mat_type != 0) {
							mat = &(pt_glb->material_parameters[mat_type - 1]);
							c11 = mat->c11; c12 = mat->c12; c44 = mat->c44;

							sxx = exxt0_crt(id) - exx0t0_crt(id);
							syy = eyyt0_crt(id) - eyy0t0_crt(id);
							szz = ezzt0_crt(id) - ezz0t0_crt(id);
							syz = eyzt0_crt(id) - eyz0t0_crt(id);
							sxz = exzt0_crt(id) - exz0t0_crt(id);
							sxy = exyt0_crt(id) - exy0t0_crt(id);

							e_sum_local = e_sum_local + \
								0.5 * c11 * (sxx * sxx + syy * syy + szz * szz) \
								+ c12 * (sxx * syy + sxx * szz + syy * szz) \
								+ 2. * c44 * (syz * syz + sxz * sxz + sxy * sxy);
						}
					}
					//get_Eelas_sum(/*exx0t0_crt, eyy0t0_crt, ezz0t0_crt, exy0t0_crt, exz0t0_crt, eyz0t0_crt, exxt0_crt, eyyt0_crt, ezzt0_crt, exyt0_crt, exzt0_crt, eyzt0_crt,*/e_sum);
				}
				e_sum = e_sum_local;

				//------Stop Criteria -> change of elastic energy is smaller than tol_elasto*e_sum for 3 consecutive steps--//
	//#pragma acc serial default(present)
	//			{
				e_elas[0] = e_elas[1]; e_elas[1] = e_elas[2]; e_elas[2] = e_elas[3]; e_elas[3] = e_sum;
				tolerance = abs(pt_glb->elasto_solver_tol * e_sum);
				//			}
				//#pragma acc wait

	//#pragma acc update host(e_elas[0:4],tolerance)
	//#pragma acc wait

				if (abs(e_elas[1] - e_elas[0]) < tolerance && \
					abs(e_elas[2] - e_elas[1]) < tolerance && \
					abs(e_elas[3] - e_elas[2]) < tolerance && \
					itrcount > 2) {
					break;
				}
			}

			if (pt_glb->if_elastostatic_film == true) {
#pragma acc parallel default(present)
				{
#pragma acc loop gang vector private(i,j,k)
					for (unsigned long id = 0; id < n; id++) {
						i = id / (ny * nz);
						j = (id - i * (ny * nz)) / nz;
						k = id - i * (ny * nz) - j * nz;

						exxt0_crt(id) = exxt0_crt(id) + exx_ext_crt(id); //+ exx0_ave_film[k];
						eyyt0_crt(id) = eyyt0_crt(id) + eyy_ext_crt(id); //+ eyy0_ave_film[k];
						ezzt0_crt(id) = ezzt0_crt(id) + ezz_ext_crt(id); //+ ezz0_ave_film[k];
						eyzt0_crt(id) = eyzt0_crt(id) + eyz_ext_crt(id); //+ eyz0_ave_film[k];
						exzt0_crt(id) = exzt0_crt(id) + exz_ext_crt(id); //+ exz0_ave_film[k];
						exyt0_crt(id) = exyt0_crt(id) + exy_ext_crt(id); //+ exy0_ave_film[k];
					}
				}
			}
		}
		//std::cout << "Elastostatic solver finish in " << itrcount << " iterations" << std::endl;

//#pragma acc parallel default(present)
//		{
//#pragma acc loop gang vector
//			for (long id = 0; id < n; id++) {
//				exxt0_crt(id) = exxt0_crt(id) + exx_ext_crt(id);
//				eyyt0_crt(id) = eyyt0_crt(id) + eyy_ext_crt(id);
//				ezzt0_crt(id) = ezzt0_crt(id) + ezz_ext_crt(id);
//				eyzt0_crt(id) = eyzt0_crt(id) + eyz_ext_crt(id);
//				exzt0_crt(id) = exzt0_crt(id) + exz_ext_crt(id);
//				exyt0_crt(id) = exyt0_crt(id) + exy_ext_crt(id);
//			}
//		}

//	else {
//#pragma acc parallel default(present)
//		{
//#pragma acc loop gang vector
//			for (long id = 0; id < n; id++) {
//				exxt0_crt(id) = exx_ext_crt(id);
//				eyyt0_crt(id) = eyy_ext_crt(id);
//				ezzt0_crt(id) = ezz_ext_crt(id);
//				eyzt0_crt(id) = eyz_ext_crt(id);
//				exzt0_crt(id) = exz_ext_crt(id);
//				exyt0_crt(id) = exy_ext_crt(id);
//			}
//		}
//	}

		if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc parallel default(present)
			{
				get_strain_static_glb();
			}
#pragma acc parallel default(present)
			{
				get_strain_static_glb_gradient();
			}
		}
	}

	//open(222, file = 'converge_statistics.dat', status = 'old', position = 'append')
	//	write(222, *) 'Electrostatic solver finish in ', itrcount!, sqrt(var / real(N))
	//	close(222)
//#pragma acc loop seq
//	for (long int i = 0; i < 1; i++) {
//		stress0homo_xx.nullify_device(); stress0homo_yy.nullify_device(); stress0homo_zz.nullify_device();
//		stress0homo_xy.nullify_device(); stress0homo_xz.nullify_device(); stress0homo_yz.nullify_device();
//
//		stress0_xx.nullify_device(); stress0_yy.nullify_device(); stress0_zz.nullify_device();
//		stress0_xy.nullify_device(); stress0_xz.nullify_device(); stress0_yz.nullify_device();
//
//		stress_xx.nullify_device(); stress_yy.nullify_device(); stress_zz.nullify_device();
//		stress_xy.nullify_device(); stress_xz.nullify_device(); stress_yz.nullify_device();
//
//		duxdx.nullify_device(); duxdy.nullify_device(); duxdz.nullify_device();
//		duydx.nullify_device(); duydy.nullify_device(); duydz.nullify_device();
//		duzdx.nullify_device(); duzdy.nullify_device(); duzdz.nullify_device();
//
//		//ux.nullify(); uy.nullify(); uz.nullify();
//
//		stress0homok_xx.nullify_device(); stress0homok_yy.nullify_device(); stress0homok_zz.nullify_device();
//		stress0homok_xy.nullify_device(); stress0homok_xz.nullify_device(); stress0homok_yz.nullify_device();
//
//		stressk_xx.nullify_device(); stressk_yy.nullify_device(); stressk_zz.nullify_device();
//		stressk_xy.nullify_device(); stressk_xz.nullify_device(); stressk_yz.nullify_device();
//
//		duxdxk.nullify_device(); duxdyk.nullify_device(); duxdzk.nullify_device();
//		duydxk.nullify_device(); duydyk.nullify_device(); duydzk.nullify_device();
//		duzdxk.nullify_device(); duzdyk.nullify_device(); duzdzk.nullify_device();
//	}
}

#pragma acc routine gang
void elastic_system::get_strain_static_1D() {
	double c11_local, c12_local;
	unsigned int mat_type;
	material* mat;

#pragma acc loop gang vector private(c11_local, c12_local, mat_type, mat)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			c11_local = mat->c11;
			c12_local = mat->c12;

			exxt0_crt(id) = exx_ext_crt(id) + exx0t0_crt(id);
			eyyt0_crt(id) = eyy_ext_crt(id) + eyy0t0_crt(id);
			ezzt0_crt(id) = -c12_local / c11_local * (exx_ext_crt(id) + eyy_ext_crt(id)) + ezz0t0_crt(id);
			eyzt0_crt(id) = eyz0t0_crt(id);
			exzt0_crt(id) = exz0t0_crt(id);
			exyt0_crt(id) = exy_ext_crt(id) + exy0t0_crt(id);
		}
		else {
			exxt0_crt(id) = 0.;
			eyyt0_crt(id) = 0.;
			ezzt0_crt(id) = 0.;
			eyzt0_crt(id) = 0.;
			exzt0_crt(id) = 0.;
			exyt0_crt(id) = 0.;
		}
	}
}

#pragma acc routine gang// nohost
void elastic_system::get_strain_static_glb() {
	double strain_in[3][3];
	double strain_out[3][3];
	unsigned int mat_type;

#pragma acc loop gang vector private(strain_in,strain_out,mat_type)
	for (long int id = 0; id < n; id++) {

		mat_type = pt_glb->material_cell(id);
		if (mat_type == 0) {
			exxt0_glb(id) = 0.;
			eyyt0_glb(id) = 0.;
			ezzt0_glb(id) = 0.;
			eyzt0_glb(id) = 0.;
			exzt0_glb(id) = 0.;
			exyt0_glb(id) = 0.;
		}
		else {
			strain_in[0][0] = exxt0_crt(id);
			strain_in[1][1] = eyyt0_crt(id);
			strain_in[2][2] = ezzt0_crt(id);
			strain_in[2][1] = eyzt0_crt(id);
			strain_in[1][2] = strain_in[2][1];
			strain_in[2][0] = exzt0_crt(id);
			strain_in[0][2] = strain_in[2][0];
			strain_in[1][0] = exyt0_crt(id);
			strain_in[0][1] = strain_in[1][0];

			pt_math->transform_matrix_crt2glb(strain_in, strain_out);

			exxt0_glb(id) = strain_out[0][0];
			eyyt0_glb(id) = strain_out[1][1];
			ezzt0_glb(id) = strain_out[2][2];
			eyzt0_glb(id) = strain_out[1][2];
			exzt0_glb(id) = strain_out[0][2];
			exyt0_glb(id) = strain_out[0][1];
		}

	}
}

#pragma acc routine gang// nohost
void elastic_system::get_strain_static_glb_gradient() {
	long i, j, k;
	unsigned int mat_type;
	long fwd, bwd;
	double scale;

#pragma acc loop gang vector private(i,j,k,fwd,bwd,scale,mat_type)
	for (long id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type == 0) {
			dexxt0dx_glb(id) = 0.; dexxt0dy_glb(id) = 0.; dexxt0dz_glb(id) = 0.;
			deyyt0dx_glb(id) = 0.; deyyt0dy_glb(id) = 0.; deyyt0dz_glb(id) = 0.;
			dezzt0dx_glb(id) = 0.; dezzt0dy_glb(id) = 0.; dezzt0dz_glb(id) = 0.;
			deyzt0dx_glb(id) = 0.; deyzt0dy_glb(id) = 0.; deyzt0dz_glb(id) = 0.;
			dexzt0dx_glb(id) = 0.; dexzt0dy_glb(id) = 0.; dexzt0dz_glb(id) = 0.;
			dexyt0dx_glb(id) = 0.; dexyt0dy_glb(id) = 0.; dexyt0dz_glb(id) = 0.;
		}
		else {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;
			//------------------spatial derivative along X--------------//
			if (stress_free_surfYZ(i, j, k) == stress_free_surfYZ(i + 1, j, k)) {
				if (stress_free_surfYZ(i, j, k) == true) {
					bwd = i; fwd = i; scale = 1.;
				}
				else {
					scale = 2.;
					if (i == 0) {
						bwd = nx - 1;
					}
					else {
						bwd = i - 1;
					}
					if (i == nx - 1) {
						fwd = 0;
					}
					else {
						fwd = i + 1;
					}
				}
			}
			else {
				if (stress_free_surfYZ(i, j, k) == true) {
					bwd = i; scale = 1.;
					if (i == nx - 1) {
						fwd = 0;
					}
					else {
						fwd = i + 1;
					}
				}
				else {
					fwd = i; scale = 1.;
					if (i == 0) {
						bwd = nx - 1;
					}
					else {
						bwd = i - 1;
					}
				}
			}
			dexxt0dx_glb(id) = (exxt0_glb(fwd, j, k) - exxt0_glb(bwd, j, k)) / dx / scale;
			deyyt0dx_glb(id) = (eyyt0_glb(fwd, j, k) - eyyt0_glb(bwd, j, k)) / dx / scale;
			dezzt0dx_glb(id) = (ezzt0_glb(fwd, j, k) - ezzt0_glb(bwd, j, k)) / dx / scale;
			deyzt0dx_glb(id) = (eyzt0_glb(fwd, j, k) - eyzt0_glb(bwd, j, k)) / dx / scale;
			dexzt0dx_glb(id) = (exzt0_glb(fwd, j, k) - exzt0_glb(bwd, j, k)) / dx / scale;
			dexyt0dx_glb(id) = (exyt0_glb(fwd, j, k) - exyt0_glb(bwd, j, k)) / dx / scale;

			//------------------spatial derivative along Y--------------//
			if (stress_free_surfXZ(i, j, k) == stress_free_surfXZ(i, j + 1, k)) {
				if (stress_free_surfXZ(i, j, k) == true) {
					bwd = j; fwd = j; scale = 1.;
				}
				else {
					scale = 2.;
					if (j == 0) {
						bwd = ny - 1;
					}
					else {
						bwd = j - 1;
					}
					if (j == ny - 1) {
						fwd = 0;
					}
					else {
						fwd = j + 1;
					}
				}
			}
			else {
				if (stress_free_surfXZ(i, j, k) == true) {
					bwd = j; scale = 1.;
					if (j == ny - 1) {
						fwd = 0;
					}
					else {
						fwd = j + 1;
					}
				}
				else {
					fwd = j; scale = 1.;
					if (j == 0) {
						bwd = ny - 1;
					}
					else {
						bwd = j - 1;
					}
				}
			}
			dexxt0dy_glb(id) = (exxt0_glb(i, fwd, k) - exxt0_glb(i, bwd, k)) / dy / scale;
			deyyt0dy_glb(id) = (eyyt0_glb(i, fwd, k) - eyyt0_glb(i, bwd, k)) / dy / scale;
			dezzt0dy_glb(id) = (ezzt0_glb(i, fwd, k) - ezzt0_glb(i, bwd, k)) / dy / scale;
			deyzt0dy_glb(id) = (eyzt0_glb(i, fwd, k) - eyzt0_glb(i, bwd, k)) / dy / scale;
			dexzt0dy_glb(id) = (exzt0_glb(i, fwd, k) - exzt0_glb(i, bwd, k)) / dy / scale;
			dexyt0dy_glb(id) = (exyt0_glb(i, fwd, k) - exyt0_glb(i, bwd, k)) / dy / scale;

			//------------------spatial derivative along Z--------------//
			if (stress_free_surfXY(i, j, k) == stress_free_surfXY(i, j, k + 1)) {
				if (stress_free_surfXY(i, j, k) == true) {
					bwd = k; fwd = k; scale = 1.;
				}
				else {
					scale = 2.;
					if (k == 0) {
						bwd = nz - 1;
					}
					else {
						bwd = k - 1;
					}
					if (k == nz - 1) {
						fwd = 0;
					}
					else {
						fwd = k + 1;
					}
				}
			}
			else {
				if (stress_free_surfXY(i, j, k) == true) {
					bwd = k; scale = 1.;
					if (k == nz - 1) {
						fwd = 0;
					}
					else {
						fwd = k + 1;
					}
				}
				else {
					fwd = k; scale = 1.;
					if (k == 0) {
						bwd = nz - 1;
					}
					else {
						bwd = k - 1;
					}
				}
			}
			dexxt0dz_glb(id) = (exxt0_glb(i, j, fwd) - exxt0_glb(i, j, bwd)) / dz / scale;
			deyyt0dz_glb(id) = (eyyt0_glb(i, j, fwd) - eyyt0_glb(i, j, bwd)) / dz / scale;
			dezzt0dz_glb(id) = (ezzt0_glb(i, j, fwd) - ezzt0_glb(i, j, bwd)) / dz / scale;
			deyzt0dz_glb(id) = (eyzt0_glb(i, j, fwd) - eyzt0_glb(i, j, bwd)) / dz / scale;
			dexzt0dz_glb(id) = (exzt0_glb(i, j, fwd) - exzt0_glb(i, j, bwd)) / dz / scale;
			dexyt0dz_glb(id) = (exyt0_glb(i, j, fwd) - exyt0_glb(i, j, bwd)) / dz / scale;
		}
	}

}
