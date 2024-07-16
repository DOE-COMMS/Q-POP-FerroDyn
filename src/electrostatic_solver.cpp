#include "ferroelectric_system.h"
#include <iostream>
//void ferroelectric_system::get_E_static() {
//	long i, j, k;
//	unsigned int itrcount;
//	unsigned int mat_type;
//	material* mat;
//	double er_inhomo;
//
//	std::complex<double> vk;
//	double e_sum_local;
//
//#pragma acc host_data use_device(px_glb.matrix, py_glb.matrix, pz_glb.matrix, pxk.matrix, pyk.matrix, pzk.matrix)
//	{
//		cufftExecD2Z(plan_electro_D2Z, px_glb.matrix, (cufftDoubleComplex*)(pxk.matrix));
//		cufftExecD2Z(plan_electro_D2Z, py_glb.matrix, (cufftDoubleComplex*)(pyk.matrix));
//		cufftExecD2Z(plan_electro_D2Z, pz_glb.matrix, (cufftDoubleComplex*)(pzk.matrix));
//		//pt_math->fourier_transform3D(px_glb, pxk);
//		//pt_math->fourier_transform3D(py_glb, pyk);
//		//pt_math->fourier_transform3D(pz_glb, pzk);
//	}
//#pragma acc wait
//
//#pragma acc parallel default(present)
//	{
//#pragma acc loop gang vector private(i,j,k)
//		for (long int id = 0; id < nx * ny * nz21; id++) {
//			i = id / (ny * nz21);
//			j = (id - i * (ny * nz21)) / nz21;
//			k = id - i * (ny * nz21) - j * nz21;
//
//			ps_terms(i, j, k) = iGq_inverse(i, j, k) / e0 * \
//				(pt_glb->kx[i] * pxk(i, j, k) + pt_glb->ky[j] * pyk(i, j, k) + pt_glb->kz[k] * pzk(i, j, k)) + \
//				im_unit * iGq_inverse(i, j, k) / e0 * rhofk(i, j, k);
//		}
//	}
//
//#pragma acc parallel default(present)
//	{
//#pragma acc loop gang vector private(i,j,k)
//		for (long int id = 0; id < nx * ny * nz21; id++) {
//			i = id / (ny * nz21);
//			j = (id - i * (ny * nz21)) / nz21;
//			k = id - i * (ny * nz21) - j * nz21;
//
//			Edxk_n(i, j, k) = -im_unit * pt_glb->kx[i] * ps_terms(i, j, k);
//			Edyk_n(i, j, k) = -im_unit * pt_glb->ky[j] * ps_terms(i, j, k);
//			Edzk_n(i, j, k) = -im_unit * pt_glb->kz[k] * ps_terms(i, j, k);
//		}
//	}
//
//#pragma acc host_data use_device(Edxk_n.matrix,Edyk_n.matrix,Edzk_n.matrix,Ex_stat.matrix,Ey_stat.matrix,Ez_stat.matrix)
//	{
//		cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edxk_n.matrix), Ex_stat.matrix);
//		cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edyk_n.matrix), Ey_stat.matrix);
//		cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edzk_n.matrix), Ez_stat.matrix);
//		//pt_math->fourier_transform3D(Edxk_n, Ex_stat); Ex_stat *= scalenn;
//		//pt_math->fourier_transform3D(Edyk_n, Ey_stat); Ey_stat *= scalenn;
//		//pt_math->fourier_transform3D(Edzk_n, Ez_stat); Ez_stat *= scalenn;
//	}
//#pragma acc wait
//
//#pragma acc parallel default(present)
//	{
//		Ex_stat *= scalenn;
//		Ey_stat *= scalenn;
//		Ez_stat *= scalenn;
//	}
//
//	e_sum_local = 0.;
//#pragma acc parallel default(present)
//	{
//#pragma acc loop gang vector reduction(+:e_sum_local)
//		for (long int id = 0; id < n; id++) {
//			e_sum_local = e_sum_local + \
//				(px_glb(id) * Ex_stat(id) + py_glb(id) * Ey_stat(id) + pz_glb(id) * Ez_stat(id)) * 0.5;
//		}
//		//	e_sum=get_Edepolar_sum(/*px_glb, py_glb, pz_glb, Ex_stat, Ey_stat, Ez_stat,*//*e_sum*/);
//	}
//	e_sum = e_sum_local;
//
//	//#pragma acc serial default(present)
//	//	{
//	e_field[0] = 0.; e_field[1] = 0.; e_field[2] = 0.; e_field[3] = 0.;
//	e_field[3] = e_sum;
//	//	}
//	//#pragma acc wait
//
//	for (itrcount = 0; itrcount < pt_glb->elec_solver_limit; itrcount++) {
//#pragma acc parallel default(present)
//		{
//#pragma acc loop gang vector private(mat_type, er_inhomo, mat)
//			for (long int id = 0; id < n; id++) {
//				mat_type = pt_glb->material_cell(id);
//				if (mat_type == 0) {
//					er_inhomo = 1. - er_homo;
//				}
//				else {
//					mat = &(pt_glb->material_parameters[mat_type - 1]);
//					er_inhomo = mat->r_permittivity - er_homo;
//				}
//				Dx_inhomo(id) = er_inhomo * Ex_stat(id);
//				Dy_inhomo(id) = er_inhomo * Ey_stat(id);
//				Dz_inhomo(id) = er_inhomo * Ez_stat(id);
//			}
//		}
//
//#pragma acc host_data use_device(Dx_inhomo.matrix, Dy_inhomo.matrix, Dz_inhomo.matrix, Dxk_inhomo.matrix, Dyk_inhomo.matrix, Dzk_inhomo.matrix)
//		{
//			cufftExecD2Z(plan_electro_D2Z, Dx_inhomo.matrix, (cufftDoubleComplex*)(Dxk_inhomo.matrix));
//			cufftExecD2Z(plan_electro_D2Z, Dy_inhomo.matrix, (cufftDoubleComplex*)(Dyk_inhomo.matrix));
//			cufftExecD2Z(plan_electro_D2Z, Dz_inhomo.matrix, (cufftDoubleComplex*)(Dzk_inhomo.matrix));
//			//pt_math->fourier_transform3D(Dx_inhomo, Dxk_inhomo);
//			//pt_math->fourier_transform3D(Dy_inhomo, Dyk_inhomo);
//			//pt_math->fourier_transform3D(Dz_inhomo, Dzk_inhomo);
//		}
//#pragma acc wait
//
//#pragma acc parallel default(present)
//		{
//#pragma acc loop gang vector private(vk,i,j,k)
//			for (long int id = 0; id < nx * ny * nz21; id++) {
//				i = id / (ny * nz21);
//				j = (id - i * (ny * nz21)) / nz21;
//				k = id - i * (ny * nz21) - j * nz21;
//
//				vk = iGq_inverse(i, j, k) * \
//					(pt_glb->kx[i] * Dxk_inhomo(i, j, k) + \
//						pt_glb->ky[j] * Dyk_inhomo(i, j, k) + \
//						pt_glb->kz[k] * Dzk_inhomo(i, j, k)) \
//					+ ps_terms(i, j, k);
//
//				Edxk_n(i, j, k) = -im_unit * pt_glb->kx[i] * vk;
//				Edyk_n(i, j, k) = -im_unit * pt_glb->ky[j] * vk;
//				Edzk_n(i, j, k) = -im_unit * pt_glb->kz[k] * vk;
//			}
//		}
//
//#pragma acc host_data use_device(Edxk_n.matrix,Edyk_n.matrix,Edzk_n.matrix,Ex_stat.matrix,Ey_stat.matrix,Ez_stat.matrix)
//		{
//			cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edxk_n.matrix), Ex_stat.matrix);
//			cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edyk_n.matrix), Ey_stat.matrix);
//			cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edzk_n.matrix), Ez_stat.matrix);
//			//pt_math->fourier_transform3D(Edxk_n, Ex_stat); Ex_stat *= scalenn;
//			//pt_math->fourier_transform3D(Edyk_n, Ey_stat); Ey_stat *= scalenn;
//			//pt_math->fourier_transform3D(Edzk_n, Ez_stat); Ez_stat *= scalenn;
//		}
//#pragma acc wait
//
//#pragma acc parallel default(present)
//		{
//			Ex_stat *= scalenn;
//			Ey_stat *= scalenn;
//			Ez_stat *= scalenn;
//		}
//
//
//		e_sum_local = 0.;
//#pragma acc parallel default(present)
//		{
//#pragma acc loop gang vector reduction(+:e_sum_local)
//			for (long int id = 0; id < n; id++) {
//				e_sum_local = e_sum_local + \
//					(px_glb(id) * Ex_stat(id) + py_glb(id) * Ey_stat(id) + pz_glb(id) * Ez_stat(id)) * 0.5;
//			}
//			//		e_sum = get_Edepolar_sum(/*px_glb, py_glb, pz_glb, Ex_stat, Ey_stat, Ez_stat,*//*e_sum*/);
//		}
//		e_sum = e_sum_local;
//
//		//------Stop Criteria -> change of depolarizing energy is smaller than tol_elec*e_sum for 3 consecutive steps--//
////#pragma acc serial default(present)
////		{
//		e_field[0] = e_field[1]; e_field[1] = e_field[2]; e_field[2] = e_field[3]; e_field[3] = e_sum;
//		tolerance = abs(pt_glb->elec_solver_tol * e_sum);
//		//		}
//		//#pragma acc wait
//
////#pragma acc update host(e_field[0:4],tolerance)
////#pragma acc wait
//
//		if (abs(e_field[1] - e_field[0]) < tolerance && \
//			abs(e_field[2] - e_field[1]) < tolerance && \
//			abs(e_field[3] - e_field[2]) < tolerance && \
//			itrcount > 2) {
//			break;
//		}
//	}
//	std::cout << "Electrostatic solver finish in " << itrcount << " iterations" << std::endl;
//	//open(222, file = 'converge_statistics.dat', status = 'old', position = 'append')
//	//	write(222, *) 'Electrostatic solver finish in ', itrcount!, sqrt(var / real(N))
//	//	close(222)
//
//	//vk.nullify();
//	//#pragma acc loop seq
//	//	for (long int i = 0; i < 1; i++) {
//	//		pxk.nullify(); pyk.nullify(); pzk.nullify();
//	//		Dxk_inhomo.nullify(); Dyk_inhomo.nullify(); Dzk_inhomo.nullify();
//	//		ps_terms.nullify();
//	//		Edxk_n.nullify(); Edyk_n.nullify(); Edzk_n.nullify();
//	//		Dx_inhomo.nullify(); Dy_inhomo.nullify(); Dz_inhomo.nullify();
//	//	}
//	//#pragma acc wait
//}

void ferroelectric_system::get_E_static/*_4RK*/() {
	long i, j, k;
	unsigned int itrcount;
	unsigned int mat_type;
	material* mat;
	double er_inhomo11, er_inhomo12, er_inhomo13;
	double er_inhomo21, er_inhomo22, er_inhomo23;
	double er_inhomo31, er_inhomo32, er_inhomo33;

	std::complex<double> vk;
	double e_sum_local;

#pragma acc host_data use_device(px_glb_store.matrix, py_glb_store.matrix, pz_glb_store.matrix, pxk.matrix, pyk.matrix, pzk.matrix)
	{
		cufftExecD2Z(plan_electro_D2Z, px_glb_store.matrix, (cufftDoubleComplex*)(pxk.matrix));
		cufftExecD2Z(plan_electro_D2Z, py_glb_store.matrix, (cufftDoubleComplex*)(pyk.matrix));
		cufftExecD2Z(plan_electro_D2Z, pz_glb_store.matrix, (cufftDoubleComplex*)(pzk.matrix));
		//pt_math->fourier_transform3D(px_glb, pxk);
		//pt_math->fourier_transform3D(py_glb, pyk);
		//pt_math->fourier_transform3D(pz_glb, pzk);
	}
#pragma acc wait

#pragma acc parallel default(present)
	{
#pragma acc loop gang vector private(i,j,k)
		for (long int id = 0; id < nx * ny * nz21; id++) {
			i = id / (ny * nz21);
			j = (id - i * (ny * nz21)) / nz21;
			k = id - i * (ny * nz21) - j * nz21;

			ps_terms(i, j, k) = iGq_inverse(i, j, k) / e0 * \
				(pt_glb->kx[i] * pxk(i, j, k) + pt_glb->ky[j] * pyk(i, j, k) + pt_glb->kz[k] * pzk(i, j, k)) + \
				im_unit * iGq_inverse(i, j, k) / e0 * rhofk(i, j, k);
		}
	}

#pragma acc parallel default(present)
	{
#pragma acc loop gang vector private(i,j,k)
		for (long int id = 0; id < nx * ny * nz21; id++) {
			i = id / (ny * nz21);
			j = (id - i * (ny * nz21)) / nz21;
			k = id - i * (ny * nz21) - j * nz21;

			Edxk_n(i, j, k) = -im_unit * pt_glb->kx[i] * ps_terms(i, j, k);
			Edyk_n(i, j, k) = -im_unit * pt_glb->ky[j] * ps_terms(i, j, k);
			Edzk_n(i, j, k) = -im_unit * pt_glb->kz[k] * ps_terms(i, j, k);
		}
	}

#pragma acc host_data use_device(Edxk_n.matrix,Edyk_n.matrix,Edzk_n.matrix,Ex_stat.matrix,Ey_stat.matrix,Ez_stat.matrix)
	{
		cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edxk_n.matrix), Ex_stat.matrix);
		cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edyk_n.matrix), Ey_stat.matrix);
		cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edzk_n.matrix), Ez_stat.matrix);
		//pt_math->fourier_transform3D(Edxk_n, Ex_stat); Ex_stat *= scalenn;
		//pt_math->fourier_transform3D(Edyk_n, Ey_stat); Ey_stat *= scalenn;
		//pt_math->fourier_transform3D(Edzk_n, Ez_stat); Ez_stat *= scalenn;
	}
#pragma acc wait

#pragma acc parallel default(present)
	{
		Ex_stat *= scalenn;
		Ey_stat *= scalenn;
		Ez_stat *= scalenn;
	}

	e_sum_local = 0.;
#pragma acc parallel default(present)
	{
#pragma acc loop gang vector reduction(+:e_sum_local)
		for (long int id = 0; id < n; id++) {
			e_sum_local = e_sum_local + \
				(px_glb_store(id) * Ex_stat(id) + py_glb_store(id) * Ey_stat(id) + pz_glb_store(id) * Ez_stat(id)) * 0.5;
		}
		//	e_sum=get_Edepolar_sum(/*px_glb, py_glb, pz_glb, Ex_stat, Ey_stat, Ez_stat,*//*e_sum*/);
	}
	e_sum = e_sum_local;

	//#pragma acc serial default(present)
	//	{
	e_field[0] = 0.; e_field[1] = 0.; e_field[2] = 0.; e_field[3] = 0.;
	e_field[3] = e_sum;
	//	}
	//#pragma acc wait

	for (itrcount = 0; itrcount < pt_glb->elec_solver_limit; itrcount++) {
#pragma acc parallel default(present)
		{
#pragma acc loop gang vector private(mat_type,\
	er_inhomo11, er_inhomo12, er_inhomo13, \
	er_inhomo21, er_inhomo22, er_inhomo23, \
	er_inhomo31, er_inhomo32, er_inhomo33, \
	mat)
			for (long int id = 0; id < n; id++) {
				mat_type = pt_glb->material_cell(id);
				if (mat_type == 0) {
					er_inhomo11 = 1. - er_homo11; er_inhomo12 = 0. - er_homo12; er_inhomo13 = 0. - er_homo13;
					er_inhomo21 = 0. - er_homo21; er_inhomo22 = 1. - er_homo22; er_inhomo23 = 0. - er_homo23;
					er_inhomo31 = 0. - er_homo31; er_inhomo32 = 0. - er_homo32; er_inhomo33 = 1. - er_homo33;
				}
				else {
					mat = &(pt_glb->material_parameters[mat_type - 1]);
					er_inhomo11 = mat->r_permittivity11 - er_homo11;
					er_inhomo12 = mat->r_permittivity12 - er_homo12;
					er_inhomo13 = mat->r_permittivity13 - er_homo13;

					er_inhomo21 = mat->r_permittivity21 - er_homo21;
					er_inhomo22 = mat->r_permittivity22 - er_homo22;
					er_inhomo23 = mat->r_permittivity23 - er_homo23;

					er_inhomo31 = mat->r_permittivity31 - er_homo31;
					er_inhomo32 = mat->r_permittivity32 - er_homo32;
					er_inhomo33 = mat->r_permittivity33 - er_homo33;
				}
				Dx_inhomo(id) = er_inhomo11 * Ex_stat(id) + er_inhomo12 * Ey_stat(id) + er_inhomo13 * Ez_stat(id);
				Dy_inhomo(id) = er_inhomo21 * Ex_stat(id) + er_inhomo22 * Ey_stat(id) + er_inhomo23 * Ez_stat(id);
				Dz_inhomo(id) = er_inhomo31 * Ex_stat(id) + er_inhomo32 * Ey_stat(id) + er_inhomo33 * Ez_stat(id);
			}
		}

#pragma acc host_data use_device(Dx_inhomo.matrix, Dy_inhomo.matrix, Dz_inhomo.matrix, Dxk_inhomo.matrix, Dyk_inhomo.matrix, Dzk_inhomo.matrix)
		{
			cufftExecD2Z(plan_electro_D2Z, Dx_inhomo.matrix, (cufftDoubleComplex*)(Dxk_inhomo.matrix));
			cufftExecD2Z(plan_electro_D2Z, Dy_inhomo.matrix, (cufftDoubleComplex*)(Dyk_inhomo.matrix));
			cufftExecD2Z(plan_electro_D2Z, Dz_inhomo.matrix, (cufftDoubleComplex*)(Dzk_inhomo.matrix));
			//pt_math->fourier_transform3D(Dx_inhomo, Dxk_inhomo);
			//pt_math->fourier_transform3D(Dy_inhomo, Dyk_inhomo);
			//pt_math->fourier_transform3D(Dz_inhomo, Dzk_inhomo);
		}
#pragma acc wait

#pragma acc parallel default(present)
		{
#pragma acc loop gang vector private(vk,i,j,k)
			for (long int id = 0; id < nx * ny * nz21; id++) {
				i = id / (ny * nz21);
				j = (id - i * (ny * nz21)) / nz21;
				k = id - i * (ny * nz21) - j * nz21;

				vk = iGq_inverse(i, j, k) * \
					(pt_glb->kx[i] * Dxk_inhomo(i, j, k) + \
						pt_glb->ky[j] * Dyk_inhomo(i, j, k) + \
						pt_glb->kz[k] * Dzk_inhomo(i, j, k)) \
					+ ps_terms(i, j, k);

				Edxk_n(i, j, k) = -im_unit * pt_glb->kx[i] * vk;
				Edyk_n(i, j, k) = -im_unit * pt_glb->ky[j] * vk;
				Edzk_n(i, j, k) = -im_unit * pt_glb->kz[k] * vk;
			}
		}

#pragma acc host_data use_device(Edxk_n.matrix,Edyk_n.matrix,Edzk_n.matrix,Ex_stat.matrix,Ey_stat.matrix,Ez_stat.matrix)
		{
			cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edxk_n.matrix), Ex_stat.matrix);
			cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edyk_n.matrix), Ey_stat.matrix);
			cufftExecZ2D(plan_electro_Z2D, (cufftDoubleComplex*)(Edzk_n.matrix), Ez_stat.matrix);
			//pt_math->fourier_transform3D(Edxk_n, Ex_stat); Ex_stat *= scalenn;
			//pt_math->fourier_transform3D(Edyk_n, Ey_stat); Ey_stat *= scalenn;
			//pt_math->fourier_transform3D(Edzk_n, Ez_stat); Ez_stat *= scalenn;
		}
#pragma acc wait

#pragma acc parallel default(present)
		{
			Ex_stat *= scalenn;
			Ey_stat *= scalenn;
			Ez_stat *= scalenn;
		}


		e_sum_local = 0.;
#pragma acc parallel default(present)
		{
#pragma acc loop gang vector reduction(+:e_sum_local)
			for (long int id = 0; id < n; id++) {
				e_sum_local = e_sum_local + \
					(px_glb_store(id) * Ex_stat(id) + py_glb_store(id) * Ey_stat(id) + pz_glb_store(id) * Ez_stat(id)) * 0.5;
			}
			//		e_sum = get_Edepolar_sum(/*px_glb, py_glb, pz_glb, Ex_stat, Ey_stat, Ez_stat,*//*e_sum*/);
		}
		e_sum = e_sum_local;

		//------Stop Criteria -> change of depolarizing energy is smaller than tol_elec*e_sum for 3 consecutive steps--//
//#pragma acc serial default(present)
//		{
		e_field[0] = e_field[1]; e_field[1] = e_field[2]; e_field[2] = e_field[3]; e_field[3] = e_sum;
		tolerance = abs(pt_glb->elec_solver_tol * e_sum);
		//		}
		//#pragma acc wait

//#pragma acc update host(e_field[0:4],tolerance)
//#pragma acc wait

		if (abs(e_field[1] - e_field[0]) < tolerance && \
			abs(e_field[2] - e_field[1]) < tolerance && \
			abs(e_field[3] - e_field[2]) < tolerance && \
			itrcount > 2) {
			break;
		}
	}
	std::cout << "Electrostatic solver finish in " << itrcount << " iterations" << std::endl;
	//open(222, file = 'converge_statistics.dat', status = 'old', position = 'append')
	//	write(222, *) 'Electrostatic solver finish in ', itrcount!, sqrt(var / real(N))
	//	close(222)

	//vk.nullify();
	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		pxk.nullify(); pyk.nullify(); pzk.nullify();
	//		Dxk_inhomo.nullify(); Dyk_inhomo.nullify(); Dzk_inhomo.nullify();
	//		ps_terms.nullify();
	//		Edxk_n.nullify(); Edyk_n.nullify(); Edzk_n.nullify();
	//		Dx_inhomo.nullify(); Dy_inhomo.nullify(); Dz_inhomo.nullify();
	//	}
	//#pragma acc wait
}

#pragma acc routine gang
void ferroelectric_system::get_E_static_1D() {
	double er_local11, er_local12, er_local13;
	double er_local21, er_local22, er_local23;
	double er_local31, er_local32, er_local33;
	double Denominator;

	double inverse_er11, inverse_er12, inverse_er13;
	double inverse_er21, inverse_er22, inverse_er23;
	double inverse_er31, inverse_er32, inverse_er33;
	unsigned int mat_type;
	material* mat;
	long x,y,z;

#pragma acc loop gang vector private(mat_type,mat,\
	er_local11,er_local12,er_local13, \
	er_local21,er_local22,er_local23, \
	er_local31,er_local32,er_local33, \
	Denominator, \
	inverse_er11, inverse_er12, inverse_er13, \
	inverse_er21, inverse_er22, inverse_er23, \
	inverse_er31, inverse_er32, inverse_er33, \
	x,y,z)
	for (long i = 0; i < n; i++) {
		mat_type = pt_glb->material_cell(i);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			er_local11 = mat->r_permittivity11; er_local12 = mat->r_permittivity12; er_local13 = mat->r_permittivity13;
			er_local21 = mat->r_permittivity21; er_local22 = mat->r_permittivity22; er_local23 = mat->r_permittivity23;
			er_local31 = mat->r_permittivity31; er_local32 = mat->r_permittivity32; er_local33 = mat->r_permittivity33;
		}
		else {
			er_local11 = 1.; er_local12 = 0.; er_local13 = 0.;
			er_local21 = 0.; er_local22 = 1.; er_local23 = 0.;
			er_local31 = 0.; er_local32 = 0.; er_local33 = 1.;
		}

		Denominator = -er_local13*er_local22*er_local31 + er_local12*er_local23*er_local31 +\
			er_local13*er_local21*er_local32 - er_local11*er_local23*er_local32 -\
			er_local12*er_local21*er_local33 + er_local11*er_local22*er_local33;

		inverse_er11 = (-er_local23*er_local32 + er_local22*er_local33) / Denominator;
		inverse_er12 = ( er_local13*er_local32 - er_local12*er_local33) / Denominator;
		inverse_er13 = (-er_local13*er_local22 + er_local12*er_local23) / Denominator;
		inverse_er21 = ( er_local23*er_local31 - er_local21*er_local33) / Denominator;
		inverse_er22 = (-er_local13*er_local31 + er_local11*er_local33) / Denominator;
		inverse_er23 = ( er_local13*er_local21 - er_local11*er_local23) / Denominator;
		inverse_er31 = (-er_local22*er_local31 + er_local21*er_local32) / Denominator;
		inverse_er32 = ( er_local12*er_local31 - er_local11*er_local32) / Denominator;
		inverse_er33 = (-er_local12*er_local21 + er_local11*er_local22) / Denominator;


		if (pt_glb->if_elec_1D_compensate == false) {
			x = i / (ny * nz);
			y = (i - x * (ny * nz)) / nz;
			z = i - x * (ny * nz) - y * nz;
			Ex_stat(i) = (0.5 * pt_glb->free_charge_cells[z] - pz_glb_store(i)) / e0 * inverse_er13;
			Ey_stat(i) = (0.5 * pt_glb->free_charge_cells[z] - pz_glb_store(i)) / e0 * inverse_er23;
			Ez_stat(i) = (0.5 * pt_glb->free_charge_cells[z] - pz_glb_store(i)) / e0 * inverse_er33;
		}
		else {
			Ex_stat(i) = 0.;
			Ey_stat(i) = 0.;
			//Ez_stat(i) = -1./(e0*er_local) * (pz_glb_store(i)-pz_t0_glb(i));
			Ez_stat(i) = 0.;
		}
	}

	//------------------ Serial method (slow)-----------------//

	////------------First cell----------//
	//mat_type = pt_glb->material_cell(0, 0, 0);
	//if (mat_type != 0) {
	//	mat = &(pt_glb->material_parameters[(mat_type)-1]);
	//	er2 = mat->r_permittivity;
	//}
	//else {
	//	er2 = 1.;
	//}

	//for (long i = 0; i < nx; i++) {
	//	for (long j = 0; j < ny; j++) {
	//		Ex_stat(i, j, 0) = 0.;
	//		Ey_stat(i, j, 0) = 0.;
	//		Ez_stat(i, j, 0) = (pt_glb->free_charge_surfaces[0] - pz_glb_store(i, j, 0) + e0 * Estat0_1D) / e0 / er2;
	//	}
	//}

	////-----------------Other cells-----------//
	//for (long k = 1; k < nz; k++) {
	//	//previous cell
	//	mat_type = pt_glb->material_cell(0, 0, k - 1);
	//	if (mat_type != 0) {
	//		mat = &(pt_glb->material_parameters[(mat_type)-1]);
	//		er1 = mat->r_permittivity;
	//	}
	//	else {
	//		er1 = 1.;
	//	}
	//	//target cell
	//	mat_type = pt_glb->material_cell(0, 0, k);
	//	if (mat_type != 0) {
	//		mat = &(pt_glb->material_parameters[(mat_type)-1]);
	//		er2 = mat->r_permittivity;
	//	}
	//	else {
	//		er2 = 1.;
	//	}
	//	//Iteration
	//	for (long i = 0; i < nx; i++) {
	//		for (long j = 0; j < ny; j++) {
	//			Ex_stat(i, j, k) = 0.;
	//			Ey_stat(i, j, k) = 0.;
	//			Ez_stat(i, j, k) = (pt_glb->free_charge_surfaces[k] + pz_glb_store(i, j, k - 1) - pz_glb_store(i, j, k) + e0 * er1 * Ez_stat(i, j, k - 1)) / e0 / er2;
	//		}
	//	}
	//}
}

//#pragma acc routine nohost
//double ferroelectric_system::get_Edepolar_sum(/*matrix3d<double>& px, matrix3d<double>& py, matrix3d<double>& pz, \
//	matrix3d<double>& Edx, matrix3d<double>& Edy, matrix3d<double>& Edz, *//*double& e_sum*/) {
//	double e_sum_local;
//#pragma acc declare device_resident(e_sum_local)
//#pragma acc loop seq 
//	for (long int i = 0; i < 1; i++) {
//		e_sum_local = 0.;
//	}
//#pragma acc wait
//
//#pragma acc loop collapse(3) independent reduction(+:e_sum_local)
//	for (long int i = 0; i < nx; i++) {
//		for (long int j = 0; j < ny; j++) {
//			for (long int k = 0; k < nz; k++) {
//				e_sum_local = e_sum_local + \
//					(px_glb(i, j, k) * Ex_stat(i, j, k) + py_glb(i, j, k) * Ey_stat(i, j, k) + pz_glb(i, j, k) * Ez_stat(i, j, k)) * 0.5;
//			}
//		}
//	}
//#pragma acc wait
//	return e_sum_local;
////#pragma acc loop seq 
////	for (long int i = 0; i < 1; i++) {
////		e_sum = e_sum_local;
////	}
//}
