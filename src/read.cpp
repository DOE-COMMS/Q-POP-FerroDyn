#include "global.h"
#include "geometry.h"

/**
 * Reads the global parameters from a file named "system_setting.in".
 * Begins reading after first 11 lines of input. Those are read by global_parameters::readgeo().
 * The file contains various parameters that define the simulation settings.
 * The function reads the parameters and sets the corresponding variables in the global_parameters class.
 * The function also builds the transformation matrix Tc2g based on the values read from the file.
 */
void global_parameters::read_global() {
	long int num_layer_unit = geometry_parameters::geo.num_layer_unit;
	std::string aline;

	std::ifstream file("system_setting.in");
	if (file.is_open()) {
		for (long int i = 0; i < 11;) {
			std::getline(file, aline);
			if (aline.length() > 0) {
				i = i + 1;
			}
		}
		file >> step;			file.ignore(1000, '\n');
		file >> dt;				file.ignore(1000, '\n');
		file >> num_materials;	file.ignore(1000, '\n');

		pt_material_layer_unit = new unsigned int[num_layer_unit];
		for (long int i = 0; i < num_layer_unit; i++) {
			file >> pt_material_layer_unit[i];
		}
		file.ignore(1000, '\n');

		file >> std::boolalpha >> if_OrientRotate;		file.ignore(1000, '\n');
		file >> Tg2c[0][0] >> Tg2c[0][1] >> Tg2c[0][2]; file.ignore(1000, '\n');
		file >> Tg2c[1][0] >> Tg2c[1][1] >> Tg2c[1][2]; file.ignore(1000, '\n');
		file >> Tg2c[2][0] >> Tg2c[2][1] >> Tg2c[2][2]; file.ignore(1000, '\n');

		file >> theta_mag >> phi_mag; file.ignore(1000, '\n');
		file >> px_init >> py_init >> pz_init; file.ignore(1000, '\n');

		file >> std::boolalpha >> if_read_Hext_nonunif;											file.ignore(1000, '\n');
		file >> Hext_stat[0] >> Hext_stat[1] >> Hext_stat[2];									file.ignore(1000, '\n');
		file >> Hext_altr[0] >> Hext_altr[1] >> Hext_altr[2];									file.ignore(1000, '\n');
		file >> H_altr_freq[0] >> H_altr_freq[1] >> H_altr_freq[2];								file.ignore(1000, '\n');
		file >> std::boolalpha >> if_prescribe_Hext >> num_prescribe_Hext;						file.ignore(1000, '\n');

		file >> Eext_stat[0] >> Eext_stat[1] >> Eext_stat[2];									file.ignore(1000, '\n');
		file >> Eext_altr[0] >> Eext_altr[1] >> Eext_altr[2];									file.ignore(1000, '\n');
		file >> E_altr_freq[0] >> E_altr_freq[1] >> E_altr_freq[2];								file.ignore(1000, '\n');
		file >> std::boolalpha >> if_prescribe_Eext >> num_prescribe_Eext;						file.ignore(1000, '\n');

		file >> std::boolalpha >> if_magnetostatics >> if_mag_1Dmodel;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_demag_factor >> demag_fac_x >> demag_fac_y >> demag_fac_z; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_elec_1Dmodel >> if_elec_1D_compensate;						file.ignore(1000, '\n');
		file >> std::boolalpha >> if_EMdynamic >> if_input_planeEM_E;;							file.ignore(1000, '\n');
		file >> planeEM_E[0] >> planeEM_E[1] >> planeEM_E[2];									file.ignore(1000, '\n');
		file >> planeEM_E_freq[0] >> planeEM_E_freq[1] >> planeEM_E_freq[2];					file.ignore(1000, '\n');
		//file >> std::boolalpha >> if_elec_1D_compensate >> if_elec_1D_dynamical_compensate		file.ignore(1000, '\n');
		file >> std::boolalpha >> if_EM_fromP >> if_EM_fromM;									file.ignore(1000, '\n');
		file >> std::boolalpha >> if_EM_backaction_on_P >> if_EM_backaction_on_M;				file.ignore(1000, '\n');
		file >> free_charge_nzi >> free_charge_nzf >> free_charge;								file.ignore(1000, '\n');
		//file >> std::boolalpha >> if_open_circuit_BC >> OC_nzi >> OC_nzf;						file.ignore(1000, '\n');
		file >> scale_elec >> elec_solver_limit >> elec_solver_tol;								file.ignore(1000, '\n');

		file >> std::boolalpha >> if_elastodynamics >> if_elasto_backaction_from_pORm >> if_elasto_on_pORm;		file.ignore(1000, '\n');
		file >> std::boolalpha >> if_gaussian_strain_pulse >> input_strain_type >> input_strain_sourcez;		file.ignore(1000, '\n');
		file >> input_strain_component >> sigma_gauss >> amplitude_gauss >> num_strain_cycles;					file.ignore(1000, '\n');
		file >> std::boolalpha >> if_elastostatic >> if_elastostatic_1D;										file.ignore(1000, '\n');
		file >> std::boolalpha >> if_read_strain_ext >> if_elastostatic_film;									file.ignore(1000, '\n');
		file >> strain_ext_nzi >> strain_ext_nzf >> exx_ext_glb >> eyy_ext_glb >> exy_ext_glb;					file.ignore(1000, '\n');
		file >> scale_elasto >> elasto_solver_limit >> elasto_solver_tol;										file.ignore(1000, '\n');

		file >> std::boolalpha >> if_1D_ABC >> if_1D_ABC_ELAST_onlytop >> if_1D_ABC_EM_onlytop >> weighting_factor_fourth_Liao; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_PEC_XY >> if_PEC_Z;										  file.ignore(1000, '\n');

		file >> std::boolalpha >> if_flexo; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_spin_pumping; file.ignore(1000, '\n');

		file >> std::boolalpha >> if_input_m >> if_input_AFMm >> if_input_Hstat;						file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_pandq >> if_input_Estat;										file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_em >> if_input_Jp;											file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_uandv >> if_input_straint0 >> if_input_eigenstraint0_crt;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_input_elastoforce;													file.ignore(1000, '\n');

		file >> std::boolalpha >> if_output_ave >> output_step_ave;			file.ignore(1000, '\n');

		file >> std::boolalpha >> if_output_m >> output_step_m >> if_output_only_magcell >> output_step_magcell;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_AFMm >> output_step_AFMm;		file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_p >> output_step_p;				file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_em >> output_step_em >> if_output_em_onecell >> output_step_em_onecell >> em_onecell_index;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_strain >> output_step_strain;	file.ignore(1000, '\n');

		file >> std::boolalpha >> if_output_Hstat >> output_step_Hstat;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_Estat >> output_step_Estat;							file.ignore(1000, '\n');
		//file >> std::boolalpha >> if_output_straint0 >> output_step_straint0;					file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_eigenstraint0_crt >> output_step_eigenstraint0_crt;	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_q >> output_step_q;									file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_uandv >> output_step_uandv;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_elastoforce >> output_step_elastoforce;				file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_emYee >> output_step_emYee;							file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_Jp >> output_step_Jp;								file.ignore(1000, '\n');
		file >> std::boolalpha >> if_output_Jishe >> output_step_Jishe;							file.ignore(1000, '\n');

		file >> std::boolalpha >> if_Jf_input;													file.ignore(1000, '\n');
		file >> Jfin_xi >> Jfin_yi >> Jfin_zi;													file.ignore(1000, '\n');
		file >> Jfin_xf >> Jfin_yf >> Jfin_zf;													file.ignore(1000, '\n');
		file >> Jf_input_type >> Jf_input_component >> Jf_input_amp >> Jf_input_freq >> num_Jf_cycles;			file.ignore(1000, '\n');
		file >> Jf_rotate_xy;																	file.ignore(1000, '\n');
		file >> std::boolalpha >> if_prescribed_m >> precess_angle >> precess_frequency;		file.ignore(1000, '\n');
	}
	file.close();
	//----------Build transformation matrix-------------//
	Tc2g[0][0] = Tg2c[0][0]; Tc2g[0][1] = Tg2c[1][0]; Tc2g[0][2] = Tg2c[2][0];
	Tc2g[1][0] = Tg2c[0][1]; Tc2g[1][1] = Tg2c[1][1]; Tc2g[1][2] = Tg2c[2][1];
	Tc2g[2][0] = Tg2c[0][2]; Tc2g[2][1] = Tg2c[1][2]; Tc2g[2][2] = Tg2c[2][2];

	read_materials();
	set_global();
}

/**
 * Reads material parameters from materials.in and initializes material objects.
 */
void global_parameters::read_materials() {

	double norm;

	material_parameters = new material[num_materials];

	std::ifstream file("materials.in");
	if (file.is_open()) {
		for (long int i = 0; i < num_materials; i++) {
			material mat;
			std::string aline;
			std::getline(file, aline);
			std::cout << "Reading material parameters of " << aline << std::endl;
			//while (aline.length() == 0)
			//	std::getline(file, aline);

			//Elastic parameters
			file >> mat.c11 >> mat.c12 >> mat.c44; 
			std::cout << "Material " << i << " - Elastic parameters: C11 = " << mat.c11 << ", C12 = " << mat.c12 << ", C44 = " << mat.c44 << std::endl;
			file.ignore(1000, '\n');
			
			file >> mat.density; 
			std::cout << "Material " << i << " - Density = " << mat.density << std::endl;
			file.ignore(1000, '\n');
			file >> mat.elast_mass_damping >> mat.elast_stiff_damping; file.ignore(1000, '\n');

			//Ferromagnetic parameters
			file >> std::boolalpha >> mat.if_FM; file.ignore(1000, '\n');
			file >> mat.Ms >> mat.gyro >> mat.FM_damping; file.ignore(1000, '\n');
			file >> mat.anisotropy_type; file.ignore(1000, '\n');
			file >> mat.K1 >> mat.K2; file.ignore(1000, '\n');
			file >> mat.Aex; file.ignore(1000, '\n');
			file >> mat.iDMI; file.ignore(1000, '\n');
			file >> mat.lamda100 >> mat.lamda111; file.ignore(1000, '\n');

			mat.B1 = -3. / 2. * mat.lamda100 * (mat.c11 - mat.c12);
			mat.B2 = -3. * mat.lamda111 * mat.c44;

			//Anti-Ferromagnetic parameters
			file >> std::boolalpha >> mat.if_AFM;													file.ignore(1000, '\n');
			file >> mat.Ms_AFM1 >> mat.Ms_AFM2;														file.ignore(1000, '\n');
			file >> mat.gyro_AFM1 >> mat.gyro_AFM2;													file.ignore(1000, '\n');
			file >> mat.damping_AFM1 >> mat.damping_AFM2;											file.ignore(1000, '\n');
			file >> mat.K1_AFM1 >> mat.K1_AFM2;														file.ignore(1000, '\n');

			file >> mat.uniaxis_AFM[0] >> mat.uniaxis_AFM[1] >> mat.uniaxis_AFM[2];					file.ignore(1000, '\n');
			norm = sqrt(mat.uniaxis_AFM[0] * mat.uniaxis_AFM[0] + mat.uniaxis_AFM[1] * mat.uniaxis_AFM[1] + mat.uniaxis_AFM[2] * mat.uniaxis_AFM[2]);
			mat.uniaxis_AFM[0] = mat.uniaxis_AFM[0] / norm;
			mat.uniaxis_AFM[1] = mat.uniaxis_AFM[1] / norm;
			mat.uniaxis_AFM[2] = mat.uniaxis_AFM[2] / norm;

			file >> mat.Aex_AFM1 >> mat.Aex_AFM2;													file.ignore(1000, '\n');
			file >> mat.J_AFM;																		file.ignore(1000, '\n');
			file >> mat.lamda100_AFM1 >> mat.lamda100_AFM2;											file.ignore(1000, '\n');
			file >> mat.lamda111_AFM1 >> mat.lamda111_AFM2;											file.ignore(1000, '\n');

			mat.B1_AFM1 = -3. / 2. * mat.lamda100_AFM1 * (mat.c11 - mat.c12);
			mat.B1_AFM2 = -3. / 2. * mat.lamda100_AFM2 * (mat.c11 - mat.c12);
			mat.B2_AFM1 = -3. * mat.lamda111_AFM1 * mat.c44;
			mat.B2_AFM2 = -3. * mat.lamda111_AFM2 * mat.c44;

			file >> std::boolalpha >> mat.if_spin_pump; file.ignore(1000, '\n');
			file >> mat.spin_hall_angle >> mat.spin_mix_cond >> mat.spin_diffus_length; file.ignore(1000, '\n');

			//Ferroelectric parameters
			file >> std::boolalpha >> mat.if_FE; file.ignore(1000, '\n');
			file >> mat.FE_mass >> mat.FE_damping; file.ignore(1000, '\n');
			file >> mat.a1; file.ignore(1000, '\n');
			file >> mat.a11 >> mat.a12; file.ignore(1000, '\n');
			file >> mat.a111 >> mat.a112 >> mat.a123; file.ignore(1000, '\n');
			file >> mat.a1111 >> mat.a1112 >> mat.a1122 >> mat.a1123; file.ignore(1000, '\n');
			file >> mat.G11; file.ignore(1000, '\n');
			file >> mat.Q11 >> mat.Q12 >> mat.Q44; file.ignore(1000, '\n');
			file >> mat.f11 >> mat.f12 >> mat.f44; file.ignore(1000, '\n');

			mat.F11 = ((mat.c11 + mat.c12) * mat.f11 - 2. * mat.c12 * mat.f12) / (mat.c11 + 2. * mat.c12) / (mat.c11 - mat.c12);
			mat.F12 = (mat.c11 * mat.f12 - mat.c12 * mat.f11) / (mat.c11 + 2. * mat.c12) / (mat.c11 - mat.c12);
			mat.F44 = mat.f44 / mat.c44;

			mat.tv1 = 2. * (mat.c11 * mat.Q11 + 2. * mat.c12 * mat.Q12);
			mat.tv2 = 2. * (mat.c11 * mat.Q12 + mat.c12 * mat.Q11 + mat.c12 * mat.Q12);
			mat.tv3 = 4. * mat.c44 * mat.Q44;
			mat.tv4 = 2. * (mat.c11 * (mat.Q11 * mat.F11 + 2. * mat.Q12 * mat.F12) + \
				2. * mat.c12 * (mat.Q11 * mat.F12 + mat.Q12 * mat.F11 + mat.Q12 * mat.F12));
			mat.tv5 = 2. * (mat.c11 * (mat.Q11 * mat.F12 + mat.Q12 * mat.F11 + mat.Q12 * mat.F12) + \
				mat.c12 * (mat.Q11 * mat.F11 + mat.Q11 * mat.F12 + mat.Q12 * mat.F11 + 3. * mat.Q12 * mat.F12));
			mat.tv6 = 2. * mat.c44 * mat.Q44 * mat.F44;
			mat.tv7 = mat.c11 * (mat.F11 * mat.F11 + 2. * mat.F12 * mat.F12) + \
				2. * mat.c12 * (mat.F12 * mat.F12 + 2. * mat.F11 * mat.F12);
			mat.tv8 = mat.c44 * mat.F44 * mat.F44;

			file >> mat.r_permittivity11 >> mat.r_permittivity12 >> mat.r_permittivity13; file.ignore(1000, '\n');
			file >> mat.r_permittivity21 >> mat.r_permittivity22 >> mat.r_permittivity23; file.ignore(1000, '\n');
			file >> mat.r_permittivity31 >> mat.r_permittivity32 >> mat.r_permittivity33; file.ignore(1000, '\n');
			file >> mat.conductivity11 >> mat.conductivity12 >> mat.conductivity13; file.ignore(1000, '\n');
			file >> mat.conductivity21 >> mat.conductivity22 >> mat.conductivity23; file.ignore(1000, '\n');
			file >> mat.conductivity31 >> mat.conductivity32 >> mat.conductivity33; file.ignore(1000, '\n');
			file >> mat.comp_n1 >> mat.omega_plasma_n1 >> mat.tao_e_n1; file.ignore(1000, '\n');
			file >> mat.comp_n2 >> mat.omega_plasma_n2 >> mat.tao_e_n2; file.ignore(1000, '\n');
			file >> mat.comp_n3 >> mat.omega_plasma_n3 >> mat.tao_e_n3; file.ignore(1000, '\n');

			// Logging the permittivity values
			std::cout << "Permittivity values:" << std::endl;
			std::cout << "r_permittivity11: " << mat.r_permittivity11 << std::endl;
			std::cout << "r_permittivity12: " << mat.r_permittivity12 << std::endl;
			std::cout << "r_permittivity13: " << mat.r_permittivity13 << std::endl;
			std::cout << "r_permittivity21: " << mat.r_permittivity21 << std::endl;
			std::cout << "r_permittivity22: " << mat.r_permittivity22 << std::endl;
			std::cout << "r_permittivity23: " << mat.r_permittivity23 << std::endl;
			std::cout << "r_permittivity31: " << mat.r_permittivity31 << std::endl;
			std::cout << "r_permittivity32: " << mat.r_permittivity32 << std::endl;
			std::cout << "r_permittivity33: " << mat.r_permittivity33 << std::endl;

			// Logging the conductivity values
			std::cout << "Conductivity values:" << std::endl;
			std::cout << "conductivity11: " << mat.conductivity11 << std::endl;
			std::cout << "conductivity12: " << mat.conductivity12 << std::endl;
			std::cout << "conductivity13: " << mat.conductivity13 << std::endl;
			std::cout << "conductivity21: " << mat.conductivity21 << std::endl;
			std::cout << "conductivity22: " << mat.conductivity22 << std::endl;
			std::cout << "conductivity23: " << mat.conductivity23 << std::endl;
			std::cout << "conductivity31: " << mat.conductivity31 << std::endl;
			std::cout << "conductivity32: " << mat.conductivity32 << std::endl;
			std::cout << "conductivity33: " << mat.conductivity33 << std::endl;


			check_material(mat);
			get_cijkl(mat);

			material_parameters[i] = mat;
		}
	}
	file.close();
}

/**
 * Reads the geometry parameters from the "system_setting.in" file and sets the corresponding member variables.
 *
 * Potential updates:
 *	1. Rework input reading to use json or xml
 * 	2. User-specified input files
 */
void geometry_parameters::readgeo() {
	std::ifstream file("system_setting.in");
	if (file.is_open()) {
		file >> nx_system >> ny_system >> nz_system; file.ignore(1000, '\n');
		file >> dx >> dy >> dz; file.ignore(1000, '\n');
		file >> std::boolalpha >> periodicX >> periodicY >> periodicZ; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_read_struct; file.ignore(1000, '\n');
		file >> nx_work >> ny_work; file.ignore(1000, '\n');
		file >> num_layer_unit; file.ignore(1000, '\n');
		file >> num_periods >> id_firstlayer >> id_lastlayer_unit; file.ignore(1000, '\n');

		pt_nz_layer_unit = new unsigned int[num_layer_unit];

		for (long int i = 0; i < num_layer_unit; i++) {
			file >> pt_nz_layer_unit[i];
		}
		file.ignore(1000, '\n');
		
		file >> std::boolalpha >> if_PML; file.ignore(1000, '\n');
		file >> std::boolalpha >> if_PML_Xs >> if_PML_Xe >> if_PML_Ys >> if_PML_Ye >> if_PML_Zs >> if_PML_Ze; file.ignore(1000, '\n');
		file >> PML_size >> sigmaMax; file.ignore(1000, '\n');
	}
	file.close();

	set_geometry();
}

void global_parameters::read_prescribe_Eext() {
	std::ifstream file("Eext.in");
	if (file.is_open()) {
		for (long int i = 0; i < num_prescribe_Eext; i++) {
			file >> prescribe_Eext_time[i] >> prescribe_Ex[i] >> prescribe_Ey[i] >> prescribe_Ez[i];
		}
	}
	file.close();
}

void global_parameters::read_prescribe_Hext() {
	std::ifstream file("Hext.in");
	if (file.is_open()) {
		for (long int i = 0; i < num_prescribe_Hext; i++) {
			file >> prescribe_Hext_time[i] >> prescribe_Hx[i] >> prescribe_Hy[i] >> prescribe_Hz[i];
		}
	}
	file.close();
}

void global_parameters::read_Hext_nonunif() {
	unsigned long x = 0;
	unsigned long y = 0;
	unsigned long z = 0;
	double Hx = 0.;
	double Hy = 0.;
	double Hz = 0.;

	std::ifstream file("Hext_nonunif.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> Hx >> Hy >> Hz;
			Hextx_nonunif(x - 1, y - 1, z - 1) = Hx;
			Hexty_nonunif(x - 1, y - 1, z - 1) = Hy;
			Hextz_nonunif(x - 1, y - 1, z - 1) = Hz;
		}
	}
	file.close();
}

void global_parameters::read_struct() {
	unsigned long x = 0;
	unsigned long y = 0;
	unsigned long z = 0;
	unsigned long id = 0;

	std::ifstream file("struct.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> id;
			material_cell(x - 1, y - 1, z - 1) = id;
		}
	}
	file.close();
}

/**
 * Logs the global parameters to a file.
 *
 * This function opens a file named "global_parameters.log" and writes the values
 * of various global parameters to it. The parameters logged include elastic,
 * ferro- and anti-ferromagnetic parameters, as well as ferroelectric parameters.
 *
 * The function does not take any parameters.
 *
 * @throws std::runtime_error if the log file cannot be opened.
 */
void global_parameters::log_global() {
	std::ofstream logFile("global_parameters.log");
	if (!logFile.is_open()) {
		std::cerr << "Unable to open log file." << std::endl;
		return;
	}

	logFile << "Global Parameters:" << std::endl;
	logFile << "nx: " << nx << std::endl;
	logFile << "ny: " << ny << std::endl;
	logFile << "nz: " << nz << std::endl;
	logFile << "n: " << n << std::endl;

	logFile << "step: " << step << std::endl;
	logFile << "dt: " << dt << std::endl;
	logFile << "time_device: " << time_device << std::endl;

	logFile << "prescribe_index_E: " << prescribe_index_E << std::endl;
	logFile << "prescribe_index_H: " << prescribe_index_H << std::endl;

	logFile << "num_materials: " << num_materials << std::endl;
	logFile << "NFE: " << NFE << std::endl;
	logFile << "NFM: " << NFM << std::endl;
	logFile << "NAFM: " << NAFM << std::endl;

	logFile << "Hext: " << Hext[0] << " " << Hext[1] << " " << Hext[2] << std::endl;
	logFile << "Hext_stat: " << Hext_stat[0] << " " << Hext_stat[1] << " " << Hext_stat[2] << std::endl;
	logFile << "Hext_altr: " << Hext_altr[0] << " " << Hext_altr[1] << " " << Hext_altr[2] << std::endl;
	logFile << "H_altr_freq: " << H_altr_freq[0] << " " << H_altr_freq[1] << " " << H_altr_freq[2] << std::endl;
	logFile << "if_read_Hext_nonunif: " << std::boolalpha << if_read_Hext_nonunif << std::endl;
	logFile << "if_prescribe_Hext: " << std::boolalpha << if_prescribe_Hext << std::endl;
	logFile << "num_prescribe_Hext: " << num_prescribe_Hext << std::endl;
	
	logFile << "Eext: " << Eext[0] << " " << Eext[1] << " " << Eext[2] << std::endl;
	logFile << "Eext_stat: " << Eext_stat[0] << " " << Eext_stat[1] << " " << Eext_stat[2] << std::endl;
	logFile << "Eext_altr: " << Eext_altr[0] << " " << Eext_altr[1] << " " << Eext_altr[2] << std::endl;
	logFile << "E_altr_freq: " << E_altr_freq[0] << " " << E_altr_freq[1] << " " << E_altr_freq[2] << std::endl;
	logFile << "if_prescribe_Eext: " << std::boolalpha << if_prescribe_Eext << std::endl;
	logFile << "num_prescribe_Eext: " << num_prescribe_Eext << std::endl;

	logFile << "if_magnetostatics: " << std::boolalpha << if_magnetostatics << std::endl;
	logFile << "if_mag_1Dmodel: " << std::boolalpha << if_mag_1Dmodel << std::endl;
	logFile << "if_elec_1Dmodel: " << std::boolalpha << if_elec_1Dmodel << std::endl;
	logFile << "if_demag_factor: " << std::boolalpha << if_demag_factor << std::endl;
	logFile << "if_elec_1D_compensate: " << std::boolalpha << if_elec_1D_compensate << std::endl;
	logFile << "if_EM_backaction_on_P: " << std::boolalpha << if_EM_backaction_on_P << std::endl;
	logFile << "if_EM_backaction_on_M: " << std::boolalpha << if_EM_backaction_on_M << std::endl;
	logFile << "if_EMdynamic: " << std::boolalpha << if_EMdynamic << std::endl;
	logFile << "if_input_planeEM_E: " << std::boolalpha << if_input_planeEM_E << std::endl;
	
	logFile << "planeEM_E: " << planeEM_E[0] << " " << planeEM_E[1] << " " << planeEM_E[2] << std::endl;
	logFile << "planeEM_E_freq: " << planeEM_E_freq[0] << " " << planeEM_E_freq[1] << " " << planeEM_E_freq[2] << std::endl;
	logFile << "if_EM_fromP: " << std::boolalpha << if_EM_fromP << std::endl;
	logFile << "if_EM_fromM: " << std::boolalpha << if_EM_fromM << std::endl;
	logFile << "if_elastodynamics: " << std::boolalpha << if_elastodynamics << std::endl;
	logFile << "if_elasto_backaction_from_pORm: " << std::boolalpha << if_elasto_backaction_from_pORm << std::endl;
	logFile << "if_elasto_on_pORm: " << std::boolalpha << if_elasto_on_pORm << std::endl;

	logFile << "if_OrientRotate: " << std::boolalpha << if_OrientRotate << std::endl;
	logFile << "Tg2c: " << Tg2c[0][0] << " " << Tg2c[0][1] << " " << Tg2c[0][2] << std::endl;
	logFile << " " << Tg2c[1][0] << " " << Tg2c[1][1] << " " << Tg2c[1][2] << std::endl;
	logFile << " " << Tg2c[2][0] << " " << Tg2c[2][1] << " " << Tg2c[2][2] << std::endl;
	logFile << "Tc2g: " << Tc2g[0][0] << " " << Tc2g[0][1] << " " << Tc2g[0][2] << std::endl;
	logFile << " " << Tc2g[1][0] << " " << Tc2g[1][1] << " " << Tc2g[1][2] << std::endl;
	logFile << " " << Tc2g[2][0] << " " << Tc2g[2][1] << " " << Tc2g[2][2] << std::endl;
	logFile << "demag_fac_x: " << demag_fac_x << std::endl;
	logFile << "demag_fac_y: " << demag_fac_y << std::endl;
	logFile << "demag_fac_z: " << demag_fac_z << std::endl;

	logFile << "if_elastostatic: " << std::boolalpha << if_elastostatic << std::endl;
	logFile << "if_elastostatic_1D: " << std::boolalpha << if_elastostatic_1D << std::endl;
	logFile << "if_elastostatic_film: " << std::boolalpha << if_elastostatic_film << std::endl;

	logFile << "free_charge_nzi: " << free_charge_nzi << std::endl;
	logFile << "free_charge_nzf: " << free_charge_nzf << std::endl;
	logFile << "free_charge: " << free_charge << std::endl;

	logFile << "if_read_strain_ext: " << std::boolalpha << if_read_strain_ext << std::endl;
	logFile << "strain_ext_nzi: " << strain_ext_nzi << std::endl;
	logFile << "strain_ext_nzf: " << strain_ext_nzf << std::endl;
	logFile << "exx_ext_glb: " << exx_ext_glb << std::endl;
	logFile << "eyy_ext_glb: " << eyy_ext_glb << std::endl;
	logFile << "exy_ext_glb: " << exy_ext_glb << std::endl;

	logFile << "scale_elec: " << scale_elec << std::endl;
	logFile << "scale_elasto: " << scale_elasto << std::endl;
	logFile << "elec_solver_limit: " << elec_solver_limit << std::endl;
	logFile << "elec_solver_tol: " << elec_solver_tol << std::endl;
	logFile << "elasto_solver_limit: " << elasto_solver_limit << std::endl;
	logFile << "elasto_solver_tol: " << elasto_solver_tol << std::endl;

	logFile << "if_gaussian_strain_pulse: " << std::boolalpha << if_gaussian_strain_pulse << std::endl;
	logFile << "input_strain_type: " << input_strain_type << std::endl;
	logFile << "input_strain_sourcez: " << input_strain_sourcez << std::endl;
	logFile << "input_strain_component: " << input_strain_component << std::endl;
	logFile << "sigma_gauss: " << sigma_gauss << std::endl;
	logFile << "amplitude_gauss: " << amplitude_gauss << std::endl;
	logFile << "num_strain_cycles: " << num_strain_cycles << std::endl;

	logFile << "if_FM_all: " << std::boolalpha << if_FM_all << std::endl;
	logFile << "if_FE_all: " << std::boolalpha << if_FE_all << std::endl;
	logFile << "if_AFM_all: " << std::boolalpha << if_AFM_all << std::endl;
	logFile << "if_periodic_allsurface: " << std::boolalpha << if_periodic_allsurface << std::endl;

	logFile << "if_1D_ABC: " << std::boolalpha << if_1D_ABC << std::endl;
	logFile << "if_1D_ABC_ELAST_onlytop: " << std::boolalpha << if_1D_ABC_ELAST_onlytop << std::endl;
	logFile << "if_1D_ABC_EM_onlytop: " << std::boolalpha << if_1D_ABC_EM_onlytop << std::endl;
	logFile << "if_PEC_XY: " << std::boolalpha << if_PEC_XY << std::endl;
	logFile << "if_PEC_Z: " << std::boolalpha << if_PEC_Z << std::endl;
	logFile << "weighting_factor_fourth_Liao: " << weighting_factor_fourth_Liao << std::endl;

	logFile << "if_PML: " << std::boolalpha << if_PML << std::endl;
	logFile << "if_PML_Xs: " << std::boolalpha << if_PML_Xs << std::endl;
	logFile << "if_PML_Xe: " << std::boolalpha << if_PML_Xe << std::endl;
	logFile << "if_PML_Ys: " << std::boolalpha << if_PML_Ys << std::endl;
	logFile << "if_PML_Ye: " << std::boolalpha << if_PML_Ye << std::endl;
	logFile << "if_PML_Zs: " << std::boolalpha << if_PML_Zs << std::endl;
	logFile << "if_PML_Ze: " << std::boolalpha << if_PML_Ze << std::endl;
	logFile << "PML_size: " << PML_size << std::endl;
	logFile << "sigmaMax: " << sigmaMax << std::endl;
	logFile << "kappaMax: " << kappaMax << std::endl;
	logFile << "PML_d: " << PML_d << std::endl;
	logFile << "PML_m: " << PML_m << std::endl;
	logFile << "eta0: " << eta0 << std::endl;
	logFile << "maxReflErr: " << maxReflErr << std::endl;
	logFile << "xS: " << xS << std::endl;
	logFile << "xE: " << xE << std::endl;
	logFile << "yS: " << yS << std::endl;
	logFile << "yE: " << yE << std::endl;
	logFile << "zS: " << zS << std::endl;
	logFile << "zE: " << zE << std::endl;

	logFile << "if_flexo: " << std::boolalpha << if_flexo << std::endl;

	logFile << "if_spin_pumping: " << std::boolalpha << if_spin_pumping << std::endl;

	logFile << "if_output_ave: " << std::boolalpha << if_output_ave << std::endl;
	logFile << "output_step_ave: " << output_step_ave << std::endl;
	
	logFile << "if_output_uandv: " << std::boolalpha << if_output_uandv << std::endl;
	logFile << "output_step_uandv: " << output_step_uandv << std::endl;
	
	logFile << "if_output_m: " << std::boolalpha << if_output_m << std::endl;
	logFile << "output_step_m: " << output_step_m << std::endl;
	logFile << "if_output_AFMm: " << std::boolalpha << if_output_AFMm << std::endl;
	logFile << "output_step_AFMm: " << output_step_AFMm << std::endl;
	logFile << "if_output_p: " << std::boolalpha << if_output_p << std::endl;
	logFile << "output_step_p: " << output_step_p << std::endl;
	logFile << "if_output_q: " << std::boolalpha << if_output_q << std::endl;
	logFile << "output_step_q: " << output_step_q << std::endl;
	logFile << "if_output_em: " << std::boolalpha << if_output_em << std::endl;
	logFile << "output_step_em: " << output_step_em << std::endl;
	logFile << "if_output_strain: " << std::boolalpha << if_output_strain << std::endl;
	logFile << "output_step_strain: " << output_step_strain << std::endl;
	logFile << "if_output_only_magcell: " << std::boolalpha << if_output_only_magcell << std::endl;
	logFile << "output_step_magcell: " << output_step_magcell << std::endl;
	logFile << "if_output_em_onecell: " << std::boolalpha << if_output_em_onecell << std::endl;
	logFile << "em_onecell_index: " << em_onecell_index << std::endl;
	logFile << "output_step_em_onecell: " << output_step_em_onecell << std::endl;
	
	logFile << "if_output_Hstat: " << std::boolalpha << if_output_Hstat << std::endl;
	logFile << "output_step_Hstat: " << output_step_Hstat << std::endl;
	logFile << "if_output_Estat: " << std::boolalpha << if_output_Estat << std::endl;
	logFile << "output_step_Estat: " << output_step_Estat << std::endl;
	logFile << "if_output_eigenstraint0_crt: " << std::boolalpha << if_output_eigenstraint0_crt << std::endl;
	logFile << "output_step_eigenstraint0_crt: " << output_step_eigenstraint0_crt << std::endl;
	
	logFile << "if_output_elastoforce: " << std::boolalpha << if_output_elastoforce << std::endl;
	logFile << "output_step_elastoforce: " << output_step_elastoforce << std::endl;
	
	logFile << "if_output_emYee: " << std::boolalpha << if_output_emYee << std::endl;
	logFile << "output_step_emYee: " << output_step_emYee << std::endl;
	logFile << "if_output_Jp: " << std::boolalpha << if_output_Jp << std::endl;
	logFile << "output_step_Jp: " << output_step_Jp << std::endl;
	logFile << "if_output_Jishe: " << std::boolalpha << if_output_Jishe << std::endl;
	logFile << "output_step_Jishe: " << output_step_Jishe << std::endl;
	
	logFile << "if_input_m: " << std::boolalpha << if_input_m << std::endl;
	logFile << "if_input_AFMm: " << std::boolalpha << if_input_AFMm << std::endl;
	logFile << "if_input_pandq: " << std::boolalpha << if_input_pandq << std::endl;
	logFile << "if_input_em: " << std::boolalpha << if_input_em << std::endl;
	logFile << "if_input_uandv: " << std::boolalpha << if_input_uandv << std::endl;
	logFile << "if_input_Hstat: " << std::boolalpha << if_input_Hstat << std::endl;
	logFile << "if_input_Estat: " << std::boolalpha << if_input_Estat << std::endl;
	logFile << "if_input_straint0: " << std::boolalpha << if_input_straint0 << std::endl;
	logFile << "if_input_eigenstraint0_crt: " << std::boolalpha << if_input_eigenstraint0_crt << std::endl;
	logFile << "if_input_Jp: " << std::boolalpha << if_input_Jp << std::endl;
	logFile << "if_input_elastoforce: " << std::boolalpha << if_input_elastoforce << std::endl;

	logFile << "if_Jf_input: " << std::boolalpha << if_Jf_input << std::endl;
	logFile << "Jf_input_type: " << Jf_input_type << std::endl;
	logFile << "Jfin_xi: " << Jfin_xi << std::endl;
	logFile << "Jfin_yi: " << Jfin_yi << std::endl;
	logFile << "Jfin_zi: " << Jfin_zi << std::endl;
	logFile << "Jfin_xf: " << Jfin_xf << std::endl;
	logFile << "Jfin_yf: " << Jfin_yf << std::endl;
	logFile << "Jfin_zf: " << Jfin_zf << std::endl;
	logFile << "Jf_input_component: " << Jf_input_component << std::endl;
	logFile << "Jf_input_amp: " << Jf_input_amp << std::endl;
	logFile << "Jf_input_freq: " << Jf_input_freq << std::endl;
	logFile << "num_Jf_cycles: " << num_Jf_cycles << std::endl;
	logFile << "Jf_rotate_xy: " << Jf_rotate_xy << std::endl;

	logFile << "if_prescribed_m: " << std::boolalpha << if_prescribed_m << std::endl;
	logFile << "precess_angle: " << precess_angle << std::endl;
	logFile << "precess_frequency: " << precess_frequency << std::endl;

	logFile << "dkx: " << dkx << std::endl;
	logFile << "dky: " << dky << std::endl;
	logFile << "dkz: " << dkz << std::endl;

	logFile.close();
}

void geometry_parameters::loggeo() {
	std::ofstream logFile("geometry_parameters.log");
	if (!logFile.is_open()) {
		std::cerr << "Unable to open geometry log file." << std::endl;
		return;
	}

	logFile << "Geometry Parameters:" << std::endl;

	logFile << "pt_nz_layer_unit: ";
	for (int i = 0; i < num_layer_unit; i++)
		logFile << pt_nz_layer_unit[i] << " ";
	logFile << std::endl;
	
	logFile << "num_periods: " << num_periods << std::endl;

	logFile << "nx_work: " << nx_work << std::endl;
	logFile << "ny_work: " << ny_work << std::endl;
	logFile << "nz_work: " << nz_work << std::endl;

	logFile << "periodicX: " << std::boolalpha << periodicX << std::endl;
	logFile << "periodicY: " << std::boolalpha << periodicY << std::endl;
	logFile << "periodicZ: " << std::boolalpha << periodicZ << std::endl;

	logFile << "num_layer_unit: " << num_layer_unit << std::endl;
	logFile << "id_lastlayer_unit: " << id_lastlayer_unit << std::endl;
	logFile << "id_firstlayer: " << id_firstlayer << std::endl;
	logFile << "id_lastlayer: " << id_lastlayer << std::endl;
	logFile << "num_layer: " << num_layer << std::endl;

	logFile << "idxi_work: " << idxi_work << std::endl;
	logFile << "idyi_work: " << idyi_work << std::endl;
	logFile << "idzi_work: " << idzi_work << std::endl;

	logFile << "idxf_work: " << idxf_work << std::endl;
	logFile << "idyf_work: " << idyf_work << std::endl;
	logFile << "idzf_work: " << idzf_work << std::endl;

	logFile << "dx: " << dx << std::endl;
	logFile << "dy: " << dy << std::endl;
	logFile << "dz: " << dz << std::endl;

	logFile << "\nPML Parameters:" << std::endl;
	logFile << "if_PML: " << std::boolalpha << if_PML << std::endl;
	logFile << "if_PML_Xs: " << if_PML_Xs << std::endl;
	logFile << "if_PML_Xe: " << if_PML_Xe << std::endl;
	logFile << "if_PML_Ys: " << if_PML_Ys << std::endl;
	logFile << "if_PML_Ye: " << if_PML_Ye << std::endl;
	logFile << "if_PML_Zs: " << if_PML_Zs << std::endl;
	logFile << "if_PML_Ze: " << if_PML_Ze << std::endl;
	logFile << "PML_size: " << PML_size << std::endl;
	logFile << "sigmaMax: " << sigmaMax << std::endl;
	logFile << "PML boundaries - xS: " << xS << ", xE: " << xE << ", yS: " << yS << ", yE: " << yE << ", zS: " << zS << ", zE: " << zE << std::endl;

	// If there are other properties or arrays that need to be logged, you can add them here in a similar fashion.

	logFile.close();
}

