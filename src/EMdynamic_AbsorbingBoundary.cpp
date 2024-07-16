#include "EMdynamic_system.h"

void EMdynamic_system::transfer_pointer() {
	if (pt_glb->if_1D_ABC == false) {
#pragma acc parallel default(present) async(8)
		{
			DEx_em_t4 = DEx_em_t3;
			DEy_em_t4 = DEy_em_t3;
			DEz_em_t4 = DEz_em_t3;
		}

#pragma acc parallel default(present) async(8)
		{
			DEx_em_t3 = DEx_em_t2;
			DEy_em_t3 = DEy_em_t2;
			DEz_em_t3 = DEz_em_t2;
		}

#pragma acc parallel default(present) async(8)
		{
			DEx_em_t2 = DEx_em_t1;
			DEy_em_t2 = DEy_em_t1;
			DEz_em_t2 = DEz_em_t1;
		}
	}

#pragma acc parallel default(present) async(8)
	{
		DEx_em_t1 = DEx_em;
		DEy_em_t1 = DEy_em;
		DEz_em_t1 = DEz_em;
	}
}

void EMdynamic_system::update_DE_Boundary_half() {
	double w;
	double coeff_bot, coeff_top;
	//#pragma acc declare device_resident(w)

	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		w = pt_glb->weighting_factor_fourth_Liao;
	//	}
	//#pragma acc wait

		//--------YZ surface -------------//
	if (pt_geo->periodicX == false) {
		//---------Y component--------//
#pragma acc parallel default(present) async(8)
		{
#pragma acc loop gang 
			for (long int j = 0; j < ny; j++) {
#pragma acc loop vector private(w) 
				for (long int k = 0; k < nz + 1; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEy_em_t1(1, j, k)) > 1.e-1) {
							DEy_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(1, j, k) - \
								(6. * w + (1. - w)) * DEy_em_t2(2, j, k) + \
								(4. * w) * DEy_em_t3(3, j, k) - \
								w * DEy_em_t4(4, j, k);
						}
						else {
							DEy_em(0, j, k) = DEy_em_t1(1, j, k);
						}

						if (abs(DEy_em_t1(nx - 1, j, k)) > 1.e-1) {
							DEy_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(nx - 1, j, k) - \
								(6. * w + (1. - w)) * DEy_em_t2(nx - 2, j, k) + \
								(4. * w) * DEy_em_t3(nx - 3, j, k) - \
								w * DEy_em_t4(nx - 4, j, k);
						}
						else {
							DEy_em(nx, j, k) = DEy_em_t1(nx - 1, j, k);
						}
					}
					else {
						DEy_em_store(0, j, k) = 0.;
						DEy_em_store(nx, j, k) = 0.;
					}
				}
			}

			//----------Z component--------//
#pragma acc loop gang 
			for (long int j = 0; j < ny + 1; j++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEz_em_t1(1, j, k)) > 1.e-1) {
							DEz_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(1, j, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(2, j, k) + \
								(4. * w) * DEz_em_t3(3, j, k) - \
								w * DEz_em_t4(4, j, k);
						}
						else {
							DEz_em(0, j, k) = DEz_em_t1(1, j, k);
						}

						if (abs(DEz_em_t1(nx - 1, j, k)) > 1.e-1) {
							DEz_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(nx - 1, j, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(nx - 2, j, k) + \
								(4. * w) * DEz_em_t3(nx - 3, j, k) - \
								w * DEz_em_t4(nx - 4, j, k);
						}
						else {
							DEz_em(nx, j, k) = DEz_em_t1(nx - 1, j, k);
						}
					}
					else {
						DEz_em_store(0, j, k) = 0.;
						DEz_em_store(nx, j, k) = 0.;
					}
				}
			}
		}
	}

	//--------XZ surface -------------//
	if (pt_geo->periodicY == false) {
		//---------X component--------//
#pragma acc parallel default(present) async(8)
		{
#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz + 1; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;
						if (abs(DEx_em_t1(i, 1, k)) > 1.e-1) {
							DEx_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, 1, k) - \
								(6. * w + (1. - w)) * DEx_em_t2(i, 2, k) + \
								(4. * w) * DEx_em_t3(i, 3, k) - \
								w * DEx_em_t4(i, 4, k);
						}
						else {
							DEx_em(i, 0, k) = DEx_em_t1(i, 1, k);
						}

						if (abs(DEx_em_t1(i, ny - 1, k)) > 1.e-1) {
							DEx_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, ny - 1, k) - \
								(6. * w + (1. - w)) * DEx_em_t2(i, ny - 2, k) + \
								(4. * w) * DEx_em_t3(i, ny - 3, k) - \
								w * DEx_em_t4(i, ny - 4, k);
						}
						else {
							DEx_em(i, ny, k) = DEx_em_t1(i, ny - 1, k);
						}
					}
					else {
						DEx_em_store(i, 0, k) = 0.;
						DEx_em_store(i, ny, k) = 0.;
					}
				}
			}
			//----------Z component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx + 1; i++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEz_em_t1(i, 1, k)) > 1.e-1) {
							DEz_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, 1, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(i, 2, k) + \
								(4. * w) * DEz_em_t3(i, 3, k) - \
								w * DEz_em_t4(i, 4, k);
						}
						else {
							DEz_em(i, 0, k) = DEz_em_t1(i, 1, k);
						}

						if (abs(DEz_em_t1(i, ny - 1, k)) > 1.e-1) {
							DEz_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, ny - 1, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(i, ny - 2, k) + \
								(4. * w) * DEz_em_t3(i, ny - 3, k) - \
								w * DEz_em_t4(i, ny - 4, k);
						}
						else {
							DEz_em(i, ny, k) = DEz_em_t1(i, ny - 1, k);
						}
					}
					else {
						DEz_em_store(i, 0, k) = 0.;
						DEz_em_store(i, ny, k) = 0.;
					}
				}
			}
		}
	}

	//--------XY surface -------------//
	if (pt_geo->periodicZ == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
#pragma acc loop vector private(w,coeff_bot,coeff_top)
				for (long int j = 0; j < ny + 1; j++) {

					if (pt_glb->if_PEC_Z == false) {
						if (pt_glb->if_1D_ABC == false) {
							w = pt_glb->weighting_factor_fourth_Liao;

							if (abs(DEx_em_t1(i, j, 1)) > 1.e-1) {
								DEx_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, 1) - \
									(6. * w + (1. - w)) * DEx_em_t2(i, j, 2) + \
									(4. * w) * DEx_em_t3(i, j, 3) - \
									w * DEx_em_t4(i, j, 4);
							}
							else {
								DEx_em(i, j, 0) = DEx_em_t1(i, j, 1);
							}

							if (abs(DEx_em_t1(i, j, nz - 1)) > 1.e-1) {
								DEx_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, nz - 1) - \
									(6. * w + (1. - w)) * DEx_em_t2(i, j, nz - 2) + \
									(4. * w) * DEx_em_t3(i, j, nz - 3) - \
									w * DEx_em_t4(i, j, nz - 4);
							}
							else {
								DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1);
							}
						}
						else {
							coeff_top = (2. * dz - c / sqrt_er_top * pt_glb->dt) / (2. * dz + c / sqrt_er_top * pt_glb->dt);
							DEx_em_store(i, j, nz) = DEx_em(i, j, nz - 1) + coeff_top * (DEx_em(i, j, nz) - DEx_em_store(i, j, nz - 1));
							if (pt_glb->if_1D_ABC_EM_onlytop == false) {
								coeff_bot = (2. * dz - c / sqrt_er_bot * pt_glb->dt) / (2. * dz + c / sqrt_er_bot * pt_glb->dt);
								DEx_em_store(i, j, 0) = DEx_em(i, j, 1) + coeff_bot * (DEx_em(i, j, 0) - DEx_em_store(i, j, 1));
							}
							else {
								DEx_em_store(i, j, 0) = 0.;
							}
						}
					}
					else {
						DEx_em_store(i, j, 0) = 0.;
						DEx_em_store(i, j, nz) = 0.;
					}

				}
			}
			//----------Y component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx + 1; i++) {
#pragma acc loop vector private(w,coeff_bot,coeff_top)
				for (long int j = 0; j < ny; j++) {

					if (pt_glb->if_PEC_Z == false) {
						if (pt_glb->if_1D_ABC == false) {
							w = pt_glb->weighting_factor_fourth_Liao;

							if (abs(DEy_em_t1(i, j, 1)) > 1.e-1) {
								DEy_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, 1) - \
									(6. * w + (1. - w)) * DEy_em_t2(i, j, 2) + \
									(4. * w) * DEy_em_t3(i, j, 3) - \
									w * DEy_em_t4(i, j, 4);
							}
							else {
								DEy_em(i, j, 0) = DEy_em_t1(i, j, 1);
							}

							if (abs(DEy_em_t1(i, j, nz - 1)) > 1.e-1) {
								DEy_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, nz - 1) - \
									(6. * w + (1. - w)) * DEy_em_t2(i, j, nz - 2) + \
									(4. * w) * DEy_em_t3(i, j, nz - 3) - \
									w * DEy_em_t4(i, j, nz - 4);
							}
							else {
								DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1);
							}
						}
						else {
							coeff_top = (2. * dz - c / sqrt_er_top * pt_glb->dt) / (2. * dz + c / sqrt_er_top * pt_glb->dt);
							DEy_em_store(i, j, nz) = DEy_em(i, j, nz - 1) + coeff_top * (DEy_em(i, j, nz) - DEy_em_store(i, j, nz - 1));
							if (pt_glb->if_1D_ABC_EM_onlytop == false) {
								coeff_bot = (2. * dz - c / sqrt_er_bot * pt_glb->dt) / (2. * dz + c / sqrt_er_bot * pt_glb->dt);
								DEy_em_store(i, j, 0) = DEy_em(i, j, 1) + coeff_bot * (DEy_em(i, j, 0) - DEy_em_store(i, j, 1));
							}
							else {
								DEy_em_store(i, j, 0) = 0.;
							}
						}
					}
					else {
						DEy_em_store(i, j, 0) = 0.;
						DEy_em_store(i, j, nz) = 0.;
					}

				}
			}
		}
	}
}

void EMdynamic_system::update_DE_Boundary_full() {
	double w;
	double coeff_bot, coeff_top;
	//#pragma acc declare device_resident(w)

	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		w = pt_glb->weighting_factor_fourth_Liao;
	//	}
	//#pragma acc wait

		//--------YZ surface -------------//
	if (pt_geo->periodicX == false) {
		//---------Y component--------//
#pragma acc parallel default(present) async(8)
		{
#pragma acc loop gang 
			for (long int j = 0; j < ny; j++) {
#pragma acc loop vector private(w) 
				for (long int k = 0; k < nz + 1; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEy_em_t1(1, j, k)) > 1.e-1) {
							DEy_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(1, j, k) - \
								(6. * w + (1. - w)) * DEy_em_t2(2, j, k) + \
								(4. * w) * DEy_em_t3(3, j, k) - \
								w * DEy_em_t4(4, j, k);
						}
						else {
							DEy_em(0, j, k) = DEy_em_t1(1, j, k);
						}

						if (abs(DEy_em_t1(nx - 1, j, k)) > 1.e-1) {
							DEy_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(nx - 1, j, k) - \
								(6. * w + (1. - w)) * DEy_em_t2(nx - 2, j, k) + \
								(4. * w) * DEy_em_t3(nx - 3, j, k) - \
								w * DEy_em_t4(nx - 4, j, k);
						}
						else {
							DEy_em(nx, j, k) = DEy_em_t1(nx - 1, j, k);
						}
					}
					else {
						DEy_em_store(0, j, k) = 0.;
						DEy_em_store(nx, j, k) = 0.;
					}
				}
			}

			//----------Z component--------//
#pragma acc loop gang 
			for (long int j = 0; j < ny + 1; j++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEz_em_t1(1, j, k)) > 1.e-1) {
							DEz_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(1, j, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(2, j, k) + \
								(4. * w) * DEz_em_t3(3, j, k) - \
								w * DEz_em_t4(4, j, k);
						}
						else {
							DEz_em(0, j, k) = DEz_em_t1(1, j, k);
						}

						if (abs(DEz_em_t1(nx - 1, j, k)) > 1.e-1) {
							DEz_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(nx - 1, j, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(nx - 2, j, k) + \
								(4. * w) * DEz_em_t3(nx - 3, j, k) - \
								w * DEz_em_t4(nx - 4, j, k);
						}
						else {
							DEz_em(nx, j, k) = DEz_em_t1(nx - 1, j, k);
						}
					}
					else {
						DEz_em_store(0, j, k) = 0.;
						DEz_em_store(nx, j, k) = 0.;
					}
				}
			}
		}
	}

	//--------XZ surface -------------//
	if (pt_geo->periodicY == false) {
		//---------X component--------//
#pragma acc parallel default(present) async(8)
		{
#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz + 1; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;
						if (abs(DEx_em_t1(i, 1, k)) > 1.e-1) {
							DEx_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, 1, k) - \
								(6. * w + (1. - w)) * DEx_em_t2(i, 2, k) + \
								(4. * w) * DEx_em_t3(i, 3, k) - \
								w * DEx_em_t4(i, 4, k);
						}
						else {
							DEx_em(i, 0, k) = DEx_em_t1(i, 1, k);
						}

						if (abs(DEx_em_t1(i, ny - 1, k)) > 1.e-1) {
							DEx_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, ny - 1, k) - \
								(6. * w + (1. - w)) * DEx_em_t2(i, ny - 2, k) + \
								(4. * w) * DEx_em_t3(i, ny - 3, k) - \
								w * DEx_em_t4(i, ny - 4, k);
						}
						else {
							DEx_em(i, ny, k) = DEx_em_t1(i, ny - 1, k);
						}
					}
					else {
						DEx_em_store(i, 0, k) = 0.;
						DEx_em_store(i, ny, k) = 0.;
					}
				}
			}
			//----------Z component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx + 1; i++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEz_em_t1(i, 1, k)) > 1.e-1) {
							DEz_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, 1, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(i, 2, k) + \
								(4. * w) * DEz_em_t3(i, 3, k) - \
								w * DEz_em_t4(i, 4, k);
						}
						else {
							DEz_em(i, 0, k) = DEz_em_t1(i, 1, k);
						}

						if (abs(DEz_em_t1(i, ny - 1, k)) > 1.e-1) {
							DEz_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, ny - 1, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(i, ny - 2, k) + \
								(4. * w) * DEz_em_t3(i, ny - 3, k) - \
								w * DEz_em_t4(i, ny - 4, k);
						}
						else {
							DEz_em(i, ny, k) = DEz_em_t1(i, ny - 1, k);
						}
					}
					else {
						DEz_em_store(i, 0, k) = 0.;
						DEz_em_store(i, ny, k) = 0.;
					}
				}
			}
		}
	}

	//--------XY surface -------------//
	if (pt_geo->periodicZ == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
#pragma acc loop vector private(w,coeff_bot,coeff_top)
				for (long int j = 0; j < ny + 1; j++) {

					if (pt_glb->if_PEC_Z == false) {
						if (pt_glb->if_1D_ABC == false) {
							w = pt_glb->weighting_factor_fourth_Liao;

							if (abs(DEx_em_t1(i, j, 1)) > 1.e-1) {
								DEx_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, 1) - \
									(6. * w + (1. - w)) * DEx_em_t2(i, j, 2) + \
									(4. * w) * DEx_em_t3(i, j, 3) - \
									w * DEx_em_t4(i, j, 4);
							}
							else {
								DEx_em(i, j, 0) = DEx_em_t1(i, j, 1);
							}

							if (abs(DEx_em_t1(i, j, nz - 1)) > 1.e-1) {
								DEx_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, nz - 1) - \
									(6. * w + (1. - w)) * DEx_em_t2(i, j, nz - 2) + \
									(4. * w) * DEx_em_t3(i, j, nz - 3) - \
									w * DEx_em_t4(i, j, nz - 4);
							}
							else {
								DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1);
							}
						}
						else {
							coeff_top = (dz - c / sqrt_er_top * pt_glb->dt) / (dz + c / sqrt_er_top * pt_glb->dt);
							DEx_em_store(i, j, nz) = DEx_em(i, j, nz - 1) + coeff_top * (DEx_em(i, j, nz) - DEx_em_store(i, j, nz - 1));
							if (pt_glb->if_1D_ABC_EM_onlytop == false) {
								coeff_bot = (dz - c / sqrt_er_bot * pt_glb->dt) / (dz + c / sqrt_er_bot * pt_glb->dt);
								DEx_em_store(i, j, 0) = DEx_em(i, j, 1) + coeff_bot * (DEx_em(i, j, 0) - DEx_em_store(i, j, 1));
							}
							else {
								DEx_em_store(i, j, 0) = 0.;
							}
						}
					}
					else {
						DEx_em_store(i, j, 0) = 0.;
						DEx_em_store(i, j, nz) = 0.;
					}

				}
			}
			//----------Y component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx + 1; i++) {
#pragma acc loop vector private(w,coeff_bot,coeff_top)
				for (long int j = 0; j < ny; j++) {

					if (pt_glb->if_PEC_Z == false) {
						if (pt_glb->if_1D_ABC == false) {
							w = pt_glb->weighting_factor_fourth_Liao;

							if (abs(DEy_em_t1(i, j, 1)) > 1.e-1) {
								DEy_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, 1) - \
									(6. * w + (1. - w)) * DEy_em_t2(i, j, 2) + \
									(4. * w) * DEy_em_t3(i, j, 3) - \
									w * DEy_em_t4(i, j, 4);
							}
							else {
								DEy_em(i, j, 0) = DEy_em_t1(i, j, 1);
							}

							if (abs(DEy_em_t1(i, j, nz - 1)) > 1.e-1) {
								DEy_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, nz - 1) - \
									(6. * w + (1. - w)) * DEy_em_t2(i, j, nz - 2) + \
									(4. * w) * DEy_em_t3(i, j, nz - 3) - \
									w * DEy_em_t4(i, j, nz - 4);
							}
							else {
								DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1);
							}
						}
						else {
							coeff_top = (dz - c / sqrt_er_top * pt_glb->dt) / (dz + c / sqrt_er_top * pt_glb->dt);
							DEy_em_store(i, j, nz) = DEy_em(i, j, nz - 1) + coeff_top * (DEy_em(i, j, nz) - DEy_em_store(i, j, nz - 1));
							if (pt_glb->if_1D_ABC_EM_onlytop == false) {
								coeff_bot = (dz - c / sqrt_er_bot * pt_glb->dt) / (dz + c / sqrt_er_bot * pt_glb->dt);
								DEy_em_store(i, j, 0) = DEy_em(i, j, 1) + coeff_bot * (DEy_em(i, j, 0) - DEy_em_store(i, j, 1));
							}
							else {
								DEy_em_store(i, j, 0) = 0.;
							}
						}
					}
					else {
						DEy_em_store(i, j, 0) = 0.;
						DEy_em_store(i, j, nz) = 0.;
					}

				}
			}
		}
	}
}

void EMdynamic_system::update_DE_Boundary() {
	double w;
	double coeff_bot, coeff_top;
	//#pragma acc declare device_resident(w)

	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		w = pt_glb->weighting_factor_fourth_Liao;
	//	}
	//#pragma acc wait

		//--------YZ surface -------------//
	if (pt_geo->periodicX == false) {
		//---------Y component--------//
#pragma acc parallel default(present) async(8)
		{
#pragma acc loop gang 
			for (long int j = 0; j < ny; j++) {
#pragma acc loop vector private(w) 
				for (long int k = 0; k < nz + 1; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEy_em_t1(1, j, k)) > 1.e-1) {
							DEy_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(1, j, k) - \
								(6. * w + (1. - w)) * DEy_em_t2(2, j, k) + \
								(4. * w) * DEy_em_t3(3, j, k) - \
								w * DEy_em_t4(4, j, k);
						}
						else {
							DEy_em(0, j, k) = DEy_em_t1(1, j, k);
						}

						if (abs(DEy_em_t1(nx - 1, j, k)) > 1.e-1) {
							DEy_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(nx - 1, j, k) - \
								(6. * w + (1. - w)) * DEy_em_t2(nx - 2, j, k) + \
								(4. * w) * DEy_em_t3(nx - 3, j, k) - \
								w * DEy_em_t4(nx - 4, j, k);
						}
						else {
							DEy_em(nx, j, k) = DEy_em_t1(nx - 1, j, k);
						}
					}
					else {
						DEy_em(0, j, k) = 0.;
						DEy_em(nx, j, k) = 0.;
					}
				}
			}

			//----------Z component--------//
#pragma acc loop gang 
			for (long int j = 0; j < ny + 1; j++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEz_em_t1(1, j, k)) > 1.e-1) {
							DEz_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(1, j, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(2, j, k) + \
								(4. * w) * DEz_em_t3(3, j, k) - \
								w * DEz_em_t4(4, j, k);
						}
						else {
							DEz_em(0, j, k) = DEz_em_t1(1, j, k);
						}

						if (abs(DEz_em_t1(nx - 1, j, k)) > 1.e-1) {
							DEz_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(nx - 1, j, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(nx - 2, j, k) + \
								(4. * w) * DEz_em_t3(nx - 3, j, k) - \
								w * DEz_em_t4(nx - 4, j, k);
						}
						else {
							DEz_em(nx, j, k) = DEz_em_t1(nx - 1, j, k);
						}
					}
					else {
						DEz_em(0, j, k) = 0.;
						DEz_em(nx, j, k) = 0.;
					}
				}
			}
		}
	}

	//--------XZ surface -------------//
	if (pt_geo->periodicY == false) {
		//---------X component--------//
#pragma acc parallel default(present) async(8)
		{
#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz + 1; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;
						if (abs(DEx_em_t1(i, 1, k)) > 1.e-1) {
							DEx_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, 1, k) - \
								(6. * w + (1. - w)) * DEx_em_t2(i, 2, k) + \
								(4. * w) * DEx_em_t3(i, 3, k) - \
								w * DEx_em_t4(i, 4, k);
						}
						else {
							DEx_em(i, 0, k) = DEx_em_t1(i, 1, k);
						}

						if (abs(DEx_em_t1(i, ny - 1, k)) > 1.e-1) {
							DEx_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, ny - 1, k) - \
								(6. * w + (1. - w)) * DEx_em_t2(i, ny - 2, k) + \
								(4. * w) * DEx_em_t3(i, ny - 3, k) - \
								w * DEx_em_t4(i, ny - 4, k);
						}
						else {
							DEx_em(i, ny, k) = DEx_em_t1(i, ny - 1, k);
						}
					}
					else {
						DEx_em(i, 0, k) = 0.;
						DEx_em(i, ny, k) = 0.;
					}
				}
			}
			//----------Z component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx + 1; i++) {
#pragma acc loop vector private(w)
				for (long int k = 0; k < nz; k++) {
					if (pt_glb->if_PEC_XY == false) {
						w = pt_glb->weighting_factor_fourth_Liao;

						if (abs(DEz_em_t1(i, 1, k)) > 1.e-1) {
							DEz_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, 1, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(i, 2, k) + \
								(4. * w) * DEz_em_t3(i, 3, k) - \
								w * DEz_em_t4(i, 4, k);
						}
						else {
							DEz_em(i, 0, k) = DEz_em_t1(i, 1, k);
						}

						if (abs(DEz_em_t1(i, ny - 1, k)) > 1.e-1) {
							DEz_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, ny - 1, k) - \
								(6. * w + (1. - w)) * DEz_em_t2(i, ny - 2, k) + \
								(4. * w) * DEz_em_t3(i, ny - 3, k) - \
								w * DEz_em_t4(i, ny - 4, k);
						}
						else {
							DEz_em(i, ny, k) = DEz_em_t1(i, ny - 1, k);
						}
					}
					else {
						DEz_em(i, 0, k) = 0.;
						DEz_em(i, ny, k) = 0.;
					}
				}
			}
		}
	}

	//--------XY surface -------------//
	if (pt_geo->periodicZ == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
#pragma acc loop vector private(w,coeff_bot,coeff_top)
				for (long int j = 0; j < ny + 1; j++) {

					if (pt_glb->if_PEC_Z == false) {
						if (pt_glb->if_1D_ABC == false) {
							w = pt_glb->weighting_factor_fourth_Liao;

							if (abs(DEx_em_t1(i, j, 1)) > 1.e-1) {
								DEx_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, 1) - \
									(6. * w + (1. - w)) * DEx_em_t2(i, j, 2) + \
									(4. * w) * DEx_em_t3(i, j, 3) - \
									w * DEx_em_t4(i, j, 4);
							}
							else {
								DEx_em(i, j, 0) = DEx_em_t1(i, j, 1);
							}

							if (abs(DEx_em_t1(i, j, nz - 1)) > 1.e-1) {
								DEx_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, nz - 1) - \
									(6. * w + (1. - w)) * DEx_em_t2(i, j, nz - 2) + \
									(4. * w) * DEx_em_t3(i, j, nz - 3) - \
									w * DEx_em_t4(i, j, nz - 4);
							}
							else {
								DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1);
							}
						}
						else {
							coeff_top = (dz - c / sqrt_er_top * pt_glb->dt) / (dz + c / sqrt_er_top * pt_glb->dt);
							DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1) + coeff_top * (DEx_em_t1(i, j, nz) - DEx_em(i, j, nz - 1));
							if (pt_glb->if_1D_ABC_EM_onlytop == false) {
								coeff_bot = (dz - c / sqrt_er_bot * pt_glb->dt) / (dz + c / sqrt_er_bot * pt_glb->dt);
								DEx_em(i, j, 0) = DEx_em_t1(i, j, 1) + coeff_bot * (DEx_em_t1(i, j, 0) - DEx_em(i, j, 1));
							}
							else {
								DEx_em(i, j, 0) = 0.;
							}
						}
					}
					else {
						DEx_em(i, j, 0) = 0.;
						DEx_em(i, j, nz) = 0.;
					}

				}
			}
			//----------Y component--------//
#pragma acc loop gang
			for (long int i = 0; i < nx + 1; i++) {
#pragma acc loop vector private(w,coeff_bot,coeff_top)
				for (long int j = 0; j < ny; j++) {

					if (pt_glb->if_PEC_Z == false) {
						if (pt_glb->if_1D_ABC == false) {
							w = pt_glb->weighting_factor_fourth_Liao;

							if (abs(DEy_em_t1(i, j, 1)) > 1.e-1) {
								DEy_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, 1) - \
									(6. * w + (1. - w)) * DEy_em_t2(i, j, 2) + \
									(4. * w) * DEy_em_t3(i, j, 3) - \
									w * DEy_em_t4(i, j, 4);
							}
							else {
								DEy_em(i, j, 0) = DEy_em_t1(i, j, 1);
							}

							if (abs(DEy_em_t1(i, j, nz - 1)) > 1.e-1) {
								DEy_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, nz - 1) - \
									(6. * w + (1. - w)) * DEy_em_t2(i, j, nz - 2) + \
									(4. * w) * DEy_em_t3(i, j, nz - 3) - \
									w * DEy_em_t4(i, j, nz - 4);
							}
							else {
								DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1);
							}
						}
						else {
							coeff_top = (dz - c / sqrt_er_top * pt_glb->dt) / (dz + c / sqrt_er_top * pt_glb->dt);
							DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1) + coeff_top * (DEy_em_t1(i, j, nz) - DEy_em(i, j, nz - 1));
							if (pt_glb->if_1D_ABC_EM_onlytop == false) {
								coeff_bot = (dz - c / sqrt_er_bot * pt_glb->dt) / (dz + c / sqrt_er_bot * pt_glb->dt);
								DEy_em(i, j, 0) = DEy_em_t1(i, j, 1) + coeff_bot * (DEy_em_t1(i, j, 0) - DEy_em(i, j, 1));
							}
							else {
								DEy_em(i, j, 0) = 0.;
							}
						}
					}
					else {
						DEy_em(i, j, 0) = 0.;
						DEy_em(i, j, nz) = 0.;
					}

				}
			}
		}
	}
}
