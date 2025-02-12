#include "global.h"
#include "geometry.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

void cleanUpLine(std::string& line, std::string& varName) {
	// Get variable name and delete it from the line
	varName = line.substr(0, line.find("="));
	varName.erase(std::remove_if(varName.begin(), varName.end(), ::isspace), varName.end());

	// Clean up line to read values
	line.erase(0, line.find("="));
	line.erase(std::remove_if(line.begin(), line.end(), [](char c) { return c == '='; }), line.end());
	line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

	auto it = line.find(';');
	if (it != std::string::npos) {
		line.erase(it);
	}

	line.erase(std::remove_if(line.begin(), line.end(), [](char c) { return c == '{' || c == '}'; }), line.end());
}

void global_parameters::read_global() {
	long int num_layer_unit = geometry_parameters::geo.num_layer_unit;

	const std::string ssFileName = "system_setting.in";

	uint32_t lineNumber = 0;

	std::ifstream ssFile(ssFileName);
	if (!ssFile.is_open()) {
		std::cerr << ssFileName << " could not be opened." << std::endl;
		std::cerr << "Check if the file exists." << std::endl;
		return;
	}

	while (!ssFile.eof()) {
		std::string line;
		std::string varName;

		std::getline(ssFile, line);
		lineNumber++;

		if (line.length() > 0) {
			if (line[0] == '#') {
				continue;
			}

			cleanUpLine(line, varName);

			std::vector<std::string> val;
			std::vector<std::string> varValues;

			std::stringstream ss(line);
			std::string token;
			while (std::getline(ss, token, ',')) {
				varValues.push_back(token);
			}

			if (varValues.size() == 0 || varValues[0] == "") {
				std::cerr << "Invalid input. Check parameter " << varName << " on line " << lineNumber << std::endl;
			}
			else {
				if (varName == "NUMBER_OF_TIMESTEPS") step = std::stoi(varValues[0]);
				else if (varName == "DT") dt = std::stod(varValues[0]);
				else if (varName == "MATERIAL_TYPES") {
					num_materials = std::stoi(varValues[0]);

					if (varValues.size() != num_materials + 1) {
						std::cerr << "Number of materials assigned does not match declared number (" << num_materials << "). Check parameter " << varName << " on line " << lineNumber << " in " << ssFileName << std::endl;
						return;
					}
					else {
						for (int i = 0; i < num_materials; i++) {
							material_names.push_back(varValues[i + 1]);
						}
					}
				}
				else if (varName == "LAYER_MATERIALS") {
					if (varValues.size() != num_layer_unit) {
						std::cerr << "Number of layers assigned does not match declared number (" << num_layer_unit << "). Check parameter " << varName << " on line " << lineNumber << " in " << ssFileName << std::endl;
						return;
					}
					else {
						pt_material_layer_unit = new unsigned int[num_layer_unit];
						for (uint64_t i = 0; i < num_layer_unit; i++) {
							pt_material_layer_unit[i] = std::stoi(varValues[i]);
						}
					}
				}

				else if (varName == "ROTATION") {
					if (varValues[0] == "TRUE") if_OrientRotate = true;
					else if_OrientRotate = false;

					if (if_OrientRotate == true) {
						double phi = std::stod(varValues[1]);
						double theta = std::stod(varValues[2]);
						double psi = std::stod(varValues[3]);

						double cphi = cos(PI * phi / 180.0);
						double sphi = sin(PI * phi / 180.0);

						double ctheta = cos(PI * theta / 180.0);
						double stheta = sin(PI * theta / 180.0);

						double cpsi = cos(PI * psi / 180.0);
						double spsi = sin(PI * psi / 180.0);

						Tg2c[0][0] = (cphi * cpsi) - (ctheta * sphi * spsi);
						Tg2c[0][1] = -(cphi * cpsi * stheta) - (cphi * spsi);
						Tg2c[0][2] = stheta * sphi;
						Tg2c[1][0] = (cpsi * sphi) + (ctheta * cphi * spsi);
						Tg2c[1][1] = (cphi * cpsi * ctheta) - (sphi * spsi);
						Tg2c[1][2] = -cphi * stheta;
						Tg2c[2][0] = stheta * spsi;
						Tg2c[2][1] = cpsi * stheta;
						Tg2c[2][2] = ctheta;
					}
					else {
						for (uint64_t i = 0; i < 3; i++) {
							for (uint64_t j = 0; j < 3; j++) {
								Tg2c[i][j] = 0.0;
							}
							Tg2c[i][i] = 1.0;
						}
					}

					Tc2g[0][0] = Tg2c[0][0]; Tc2g[0][1] = Tg2c[1][0]; Tc2g[0][2] = Tg2c[2][0];
					Tc2g[1][0] = Tg2c[0][1]; Tc2g[1][1] = Tg2c[1][1]; Tc2g[1][2] = Tg2c[2][1];
					Tc2g[2][0] = Tg2c[0][2]; Tc2g[2][1] = Tg2c[1][2]; Tc2g[2][2] = Tg2c[2][2];
				}

				else if (varName == "THETA_MAG") theta_mag = std::stod(varValues[0]);
				else if (varName == "PHI_MAG") phi_mag = std::stod(varValues[0]);
				else if (varName == "INITIAL_P") {
					px_init = std::stod(varValues[0]);
					py_init = std::stod(varValues[1]);
					pz_init = std::stod(varValues[2]);
				}

				else if (varName == "READ_H_EXT_NONUNIFORM") {
					if (varValues[0] == "TRUE") if_read_Hext_nonunif = true;
					else if_read_Hext_nonunif = false;
				}
				else if (varName == "H_EXT_STATIC" && if_read_Hext_nonunif) {
					Hext_stat[0] = std::stod(varValues[0]);
					Hext_stat[1] = std::stod(varValues[1]);
					Hext_stat[2] = std::stod(varValues[2]);
				}
				else if (varName == "H_EXT_ALTERNATING" && if_read_Hext_nonunif) {
					Hext_altr[0] = std::stod(varValues[0]);
					Hext_altr[1] = std::stod(varValues[1]);
					Hext_altr[2] = std::stod(varValues[2]);
				}
				else if (varName == "H_EXT_FREQ" && if_read_Hext_nonunif) {
					H_altr_freq[0] = std::stod(varValues[0]);
					H_altr_freq[1] = std::stod(varValues[1]);
					H_altr_freq[2] = std::stod(varValues[2]);
				}
				else if (varName == "READ_H_EXT_PRESCRIBED") {
					if (varValues[0] == "TRUE") if_prescribe_Hext = true;
					else if_prescribe_Hext = false;

					num_prescribe_Hext = std::stoi(varValues[1]);
				}

				else if (varName == "E_EXT_STATIC") {
					Eext_stat[0] = std::stod(varValues[0]);
					Eext_stat[1] = std::stod(varValues[1]);
					Eext_stat[2] = std::stod(varValues[2]);
				}
				else if (varName == "E_EXT_ALTERNATING") {
					Eext_altr[0] = std::stod(varValues[0]);
					Eext_altr[1] = std::stod(varValues[1]);
					Eext_altr[2] = std::stod(varValues[2]);
				}
				else if (varName == "E_EXTFREQ") {
					E_altr_freq[0] = std::stod(varValues[0]);
					E_altr_freq[1] = std::stod(varValues[1]);
					E_altr_freq[2] = std::stod(varValues[2]);
				}
				else if (varName == "READ_E_EXT_PRESCRIBED") {
					if (varValues[0] == "TRUE") if_prescribe_Eext = true;
					else if_prescribe_Eext = false;

					num_prescribe_Eext = std::stoi(varValues[1]);
				}

				else if (varName == "MAGNETOSTATICS") {
					if (varValues[0] == "TRUE") if_magnetostatics = true;
					else if_magnetostatics = false;
				}
				else if (varName == "DEMAG_FIELD_1D") {
					if (varValues[0] == "TRUE") if_mag_1Dmodel = true;
					else if_mag_1Dmodel = false;
				}
				else if (varName == "DEMAG_FACTOR") {
					if (varValues[0] == "TRUE") if_demag_factor = true;
					else if_demag_factor = false;

					demag_fac_x = std::stod(varValues[1]);
					demag_fac_y = std::stod(varValues[2]);
					demag_fac_z = std::stod(varValues[3]);
				}

				else if (varName == "DEPOLARIZING_FIELD_1D") {
					if (varValues[0] == "TRUE") if_elec_1Dmodel = true;
					else if_elec_1Dmodel = false;
				}
				else if (varName == "ED_ZERO_FIELD") {
					if (varValues[0] == "TRUE") if_elec_1D_compensate = true;
					else if_elec_1D_compensate = false;
				}
				else if (varName == "EM_DYNAMICS") {
					if (varValues[0] == "TRUE") if_EMdynamic = true;
					else if_EMdynamic = false;
				}
				else if (varName == "INJECT_EM_PLANEWAVE") {
					if (varValues[0] == "TRUE") if_input_planeEM_E = true;
					else if_input_planeEM_E = false;
				}
				else if (varName == "INJECT_EM_PLANEWAVE_AMPL" && if_input_planeEM_E) {
					planeEM_E[0] = std::stod(varValues[0]);
					planeEM_E[1] = std::stod(varValues[1]);
					planeEM_E[2] = std::stod(varValues[2]);
				}
				else if (varName == "INJECT_EM_PLANEWAVE_FREQ" && if_input_planeEM_E) {
					planeEM_E_freq[0] = std::stod(varValues[0]);
					planeEM_E_freq[1] = std::stod(varValues[1]);
					planeEM_E_freq[2] = std::stod(varValues[2]);
				}
				else if (varName == "EM_WAVE_FROM_P") {
					if (varValues[0] == "TRUE") if_EM_fromP = true;
					else if_EM_fromP = false;
				}
				else if (varName == "EM_WAVE_FROM_M") {
					if (varValues[0] == "TRUE") if_EM_fromM = true;
					else if_EM_fromM = false;
				}
				else if (varName == "EM_BACKACTION_P") {
					if (varValues[0] == "TRUE") if_EM_backaction_on_P = true;
					else if_EM_backaction_on_P = false;
				}
				else if (varName == "EM_BACKACTION_M") {
					if (varValues[0] == "TRUE") if_EM_backaction_on_M = true;
					else if_EM_backaction_on_M = false;
				}

				else if (varName == "SURFACE_CHARGE") {
					free_charge_nzi = std::stoi(varValues[0]);
					free_charge_nzf = std::stoi(varValues[1]);
					free_charge = std::stod(varValues[2]);
				}

				else if (varName == "ELECTROSTATIC_SOLVER_SCALING") {
					scale_elec = std::stod(varValues[0]);
				}
				else if (varName == "ELECTROSTATIC_SOLVER_MAX_ITERATIONS") {
					elec_solver_limit = std::stod(varValues[0]);
				}
				else if (varName == "ELECTROSTATIC_SOLVER_TOLERANCE") {
					elec_solver_tol = std::stod(varValues[0]);
				}

				else if (varName == "ELASTODYNAMICS") {
					if (varValues[0] == "TRUE") if_elastodynamics = true;
					else if_elastodynamics = false;
				}
				else if (varName == "ELASTODYNAMICS_WAVE_ON_PM") {
					if (varValues[0] == "TRUE") if_elasto_on_pORm = true;
					else if_elasto_on_pORm = false;
				}
				else if (varName == "ELASTODYNAMICS_BACKACTION_ON_PM") {
					if (varValues[0] == "TRUE") if_elasto_backaction_from_pORm = true;
					else if_elasto_backaction_from_pORm = false;
				}

				else if (varName == "GAUSSIAN_STRAIN_PULSE") {
					if (varValues[0] == "TRUE") if_gaussian_strain_pulse = true;
					else if_gaussian_strain_pulse = false;

					input_strain_type = std::stoi(varValues[1]);
					input_strain_sourcez = std::stoi(varValues[2]);

					if (varValues[3] == "x" || varValues[3] == "X") input_strain_component = 'x';
					else if (varValues[3] == "y" || varValues[3] == "Y") input_strain_component = 'y';
					else if (varValues[3] == "z" || varValues[3] == "Z") input_strain_component = 'z';
					else if (if_gaussian_strain_pulse == true) std::cerr << "input_strain_component must be x, y, or z" << std::endl;

					sigma_gauss = std::stod(varValues[4]);
					amplitude_gauss = std::stod(varValues[5]);
					num_strain_cycles = std::stoi(varValues[6]);
				}

				else if (varName == "ELASTOSTATIC") {
					if (varValues[0] == "TRUE") if_elastostatic = true;
					else if_elastostatic = false;

					if (varValues[1] == "TRUE") if_elastostatic_1D = true;
					else if_elastostatic_1D = false;
				}

				else if (varName == "READ_STRAIN_EXT") {
					if (varValues[0] == "TRUE") if_read_strain_ext = true;
					else if_read_strain_ext = false;

					if (varValues[1] == "TRUE") if_elastostatic_film = true;
					else if_elastostatic_film = false;

					strain_ext_nzi = std::stoi(varValues[2]);
					strain_ext_nzf = std::stoi(varValues[3]);
					exx_ext_glb = std::stod(varValues[4]);
					eyy_ext_glb = std::stod(varValues[5]);
					exy_ext_glb = std::stod(varValues[6]);
				}
				else if (varName == "ELASTOSTATIC_SOLVER_SCALING") {
					scale_elasto = std::stod(varValues[0]);
				}
				else if (varName == "ELASTOSTATIC_SOLVER_MAX_ITERATIONS") {
					elasto_solver_limit = std::stod(varValues[0]);
				}
				else if (varName == "ELASTOSTATIC_SOLVER_TOLERANCE") {
					elasto_solver_tol = std::stod(varValues[0]);
				}

				else if (varName == "1D_ABC") {
					if (varValues[0] == "TRUE") if_1D_ABC = true;
					else if_1D_ABC = false;

					if (varValues[1] == "TRUE") if_1D_ABC_ELAST_onlytop = true;
					else if_1D_ABC_ELAST_onlytop = false;

					if (varValues[2] == "TRUE") if_1D_ABC_EM_onlytop = true;
					else if_1D_ABC_EM_onlytop = false;

					weighting_factor_fourth_Liao = std::stod(varValues[3]);
				}

				else if (varName == "PEC") {
					if (varValues[0] == "TRUE") if_PEC_XY = true;
					else if_PEC_XY = false;

					if (varValues[1] == "TRUE") if_PEC_Z = true;
					else if_PEC_Z = false;
				}

				else if (varName == "FLEXOELECTRICITY") {
					if (varValues[0] == "TRUE") if_flexo = true;
					else if_flexo = false;
				}

				else if (varName == "SPIN_PUMPING") {
					if (varValues[0] == "TRUE") if_spin_pumping = true;
					else if_spin_pumping = false;
				}

				else if (varName == "INPUT_M") {
					if (varValues[0] == "TRUE") if_input_m = true;
					else if_input_m = false;
				}
				else if (varName == "INPUT_AFM") {
					if (varValues[0] == "TRUE") if_input_AFMm = true;
					else if_input_AFMm = false;
				}
				else if (varName == "INPUT_HSTAT") {
					if (varValues[0] == "TRUE") if_input_Hstat = true;
					else if_input_Hstat = false;
				}
				else if (varName == "INPUT_P_Q") {
					if (varValues[0] == "TRUE") if_input_pandq = true;
					else if_input_pandq = false;
				}
				else if (varName == "INPUT_ESTAT") {
					if (varValues[0] == "TRUE") if_input_Estat = true;
					else if_input_Estat = false;
				}
				else if (varName == "INPUT_EM") {
					if (varValues[0] == "TRUE") if_input_em = true;
					else if_input_em = false;
				}
				else if (varName == "INPUT_JP") {
					if (varValues[0] == "TRUE") if_input_Jp = true;
					else if_input_Jp = false;
				}
				else if (varName == "INPUT_U_V") {
					if (varValues[0] == "TRUE") if_input_uandv = true;
					else if_input_uandv = false;
				}
				else if (varName == "INPUT_INITIAL_STRAIN") {
					if (varValues[0] == "TRUE") if_input_straint0 = true;
					else if_input_straint0 = false;
				}
				else if (varName == "INPUT_EIGENSTRAIN") {
					if (varValues[0] == "TRUE") if_input_eigenstraint0_crt = true;
					else if_input_eigenstraint0_crt = false;
				}
				else if (varName == "INPUT_ELASTOFORCE") {
					if (varValues[0] == "TRUE") if_input_elastoforce = true;
					else if_input_elastoforce = false;
				}
				else if (varName == "OUTPUT_AVERAGE_P_M") {
					if (varValues[0] == "TRUE") if_output_ave = true;
					else if_output_ave = false;

					output_step_ave = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_M") {
					if (varValues[0] == "TRUE") if_output_m = true;
					else if_output_m = false;

					output_step_m = std::stoi(varValues[1]);

					if (varValues[2] == "TRUE") if_output_only_magcell = true;
					else if_output_only_magcell = false;

					output_step_magcell = std::stoi(varValues[3]);
				}
				else if (varName == "OUTPUT_AFM_M") {
					if (varValues[0] == "TRUE") if_output_AFMm = true;
					else if_output_AFMm = false;

					output_step_AFMm = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_P") {
					if (varValues[0] == "TRUE") if_output_p = true;
					else if_output_p = false;

					output_step_p = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_EM") {
					if (varValues[0] == "TRUE") if_output_em = true;
					else if_output_em = false;

					output_step_em = std::stoi(varValues[1]);

					if (varValues[2] == "TRUE") if_output_em_onecell = true;
					else if_output_em_onecell = false;

					output_step_em_onecell = std::stoi(varValues[3]);
					em_onecell_index = std::stoi(varValues[4]);
				}
				else if (varName == "OUTPUT_STRAIN") {
					if (varValues[0] == "TRUE") if_output_strain = true;
					else if_output_strain = false;

					output_step_strain = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_HSTAT") {
					if (varValues[0] == "TRUE") if_output_Hstat = true;
					else if_output_Hstat = false;

					output_step_Hstat = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_ESTAT") {
					if (varValues[0] == "TRUE") if_output_Estat = true;
					else if_output_Estat = false;

					output_step_Estat = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_EIGENSTRAIN") {
					if (varValues[0] == "TRUE") if_output_eigenstraint0_crt = true;
					else if_output_eigenstraint0_crt = false;

					output_step_eigenstraint0_crt = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_Q") {
					if (varValues[0] == "TRUE") if_output_q = true;
					else if_output_q = false;

					output_step_q = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_U_V") {
					if (varValues[0] == "TRUE") if_output_uandv = true;
					else if_output_uandv = false;

					output_step_uandv = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_ELASTOFORCE") {
					if (varValues[0] == "TRUE") if_output_elastoforce = true;
					else if_output_elastoforce = false;

					output_step_elastoforce = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_EM_YEE") {
					if (varValues[0] == "TRUE") if_output_emYee = true;
					else if_output_emYee = false;

					output_step_emYee = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_JP") {
					if (varValues[0] == "TRUE") if_output_Jp = true;
					else if_output_Jp = false;

					output_step_Jp = std::stoi(varValues[1]);
				}
				else if (varName == "OUTPUT_JISHE") {
					if (varValues[0] == "TRUE") if_output_Jishe = true;
					else if_output_Jishe = false;

					output_step_Jishe = std::stoi(varValues[1]);
				}

				else if (varName == "IF_JF_INPUT") {
					if (varValues[0] == "TRUE") if_Jf_input = true;
					else if_Jf_input = false;
				}
				else if (varName == "JF_POSITION" && if_Jf_input) {
					if (varValues.size() == 6) {
						Jfin_xi = std::stod(varValues[0]);
						Jfin_yi = std::stod(varValues[1]);
						Jfin_zi = std::stod(varValues[2]);
						Jfin_xf = std::stod(varValues[3]);
						Jfin_yf = std::stod(varValues[4]);
						Jfin_zf = std::stod(varValues[5]);
					}
				}
				else if (varName == "JF_INPUT_TYPE" && if_Jf_input) {
					if (varValues.size() == 6) {
						Jf_input_type = std::stoi(varValues[0]);

						if (varValues[1] == "x" || varValues[1] == "X") Jf_input_component = 'x';
						else if (varValues[1] == "y" || varValues[1] == "Y") Jf_input_component = 'y';
						else if (varValues[1] == "z" || varValues[1] == "Z") Jf_input_component = 'z';
						else if (if_Jf_input == true) std::cerr << "Jf_input_component must be x, y, or z" << std::endl;

						Jf_input_amp = std::stod(varValues[2]);
						Jf_input_freq = std::stod(varValues[3]);
						num_Jf_cycles = std::stoi(varValues[4]);
						Jf_rotate_xy = std::stoi(varValues[5]);
					}
				}
				else if (varName == "IF_PRESCRIBED_M" && if_Jf_input) {
					if (varValues[0] == "TRUE") if_prescribed_m = true;
					else if_prescribed_m = false;

					precess_angle = std::stod(varValues[1]);
					precess_frequency = std::stod(varValues[2]);
				}
			}
		}
	}

	ssFile.close();

	std::cout << "Global parameters read successfully from " << ssFileName << std::endl;

	read_materials();
	set_global();
}

/**
 * Reads material parameters from materials.in and initializes material objects.
 */
void global_parameters::read_materials() {
	material_parameters = new material[num_materials];

	for (long int i = 0; i < num_materials; i++) {
		material mat;

		const std::string ssFileName = material_names.at(i) + ".in";
		std::ifstream file(ssFileName);

		int lineNumber = 0;

		std::ifstream ssFile(ssFileName);
		if (!ssFile.is_open()) {
			std::cerr << ssFileName << " could not be opened." << std::endl;
			std::cerr << "Check if the file exists." << std::endl;
			return;
		}

		while (!ssFile.eof()) {
			std::string line;
			std::string varName;

			std::getline(ssFile, line);
			lineNumber++;

			if (line.length() > 0) {
				if (line[0] == '#') {
					continue;
				}

				cleanUpLine(line, varName);

				std::vector<std::string> val;
				std::vector<std::string> varValues;

				std::stringstream ss(line);
				std::string token;
				while (std::getline(ss, token, ',')) {
					varValues.push_back(token);
				}

				if (varValues.size() == 0 || varValues[0] == "") {
					std::cerr << "Invalid input. Check parameter " << varName << " on line " << lineNumber << std::endl;
				}
				else {
					if (varName == "MATERIAL" && varValues[0] != material_names.at(i)) {
						std::cerr << "Wrong material file read for material " << material_names.at(i) << std::endl;
						std::cerr << "Check material name and parameters in " << ssFileName << std::endl;
						return;
					}
					else if (varName == "C11") mat.c11 = std::stod(varValues[0]);
					else if (varName == "C12") mat.c12 = std::stod(varValues[0]);
					else if (varName == "C44") mat.c44 = std::stod(varValues[0]);
					else if (varName == "DENSITY") mat.density = std::stod(varValues[0]);
					else if (varName == "ELASTODYNAMIC_MASS") mat.elast_mass_damping = std::stod(varValues[0]);
					else if (varName == "ELASTODYNAMIC_DAMPING") mat.elast_stiff_damping = std::stod(varValues[0]);

					else if (varName == "IS_MAGNETIC") {
						if (varValues[0] == "TRUE") mat.if_FM = true;
						else mat.if_FM = false;
					}
					else if (varName == "FM_SATURATION_MAGNETIZATION") mat.Ms = std::stod(varValues[0]);
					else if (varName == "FM_GYROMAGNETIC_RATIO") mat.gyro = std::stod(varValues[0]);
					else if (varName == "FM_DAMPING") mat.FM_damping = std::stod(varValues[0]);
					else if (varName == "FM_MAG_ANISO_TYPE") {
						if (varValues[0] == "0") mat.anisotropy_type = 0;
						else if (varValues[0] == "1") mat.anisotropy_type = 1;
						else if (varValues[0] == "2") mat.anisotropy_type = 2;
						else {
							std::cerr << "Magnetocrystalline anisotropy parameter ("<< varValues[0] <<") must be 0, 1, or 2. Defaulting to 0.\n";
							mat.anisotropy_type = 0;
						}
					}
					else if (varName == "FM_K1") mat.K1 = std::stod(varValues[0]);
					else if (varName == "FM_K2") mat.K2 = std::stod(varValues[0]);
					else if (varName == "FM_EXCHANGE_ANISO") mat.Aex = std::stod(varValues[0]);
					else if (varName == "FM_DMI_STRENGTH") mat.iDMI = std::stod(varValues[0]);
					else if (varName == "FM_LAMBDA_100") mat.lamda100 = std::stod(varValues[0]);
					else if (varName == "FM_LAMBDA_111") mat.lamda111 = std::stod(varValues[0]);

					else if (varName == "IS_ANTIFERROMAGNETIC") {
						if (varValues[0] == "TRUE") mat.if_AFM = true;
						else mat.if_AFM = false;
					}
					else if (varName == "AFM_SATURATION_MAGNETIZATION") {
						mat.Ms_AFM1 = std::stod(varValues[0]);
						mat.Ms_AFM2 = std::stod(varValues[1]);
					}
					else if (varName == "AFM_GYROMAGNETIC_RATIO") {
						mat.gyro_AFM1 = std::stod(varValues[0]);
						mat.gyro_AFM2 = std::stod(varValues[1]);
					}
					else if (varName == "AFM_DAMPING") {
						mat.damping_AFM1 = std::stod(varValues[0]);
						mat.damping_AFM2 = std::stod(varValues[1]);
					}
					else if (varName == "AFM_K1") {
						mat.K1_AFM1 = std::stod(varValues[0]);
						mat.K1_AFM2 = std::stod(varValues[1]);
					}
					else if (varName == "AFM_UNIAXIAL_ANISO") {
						mat.uniaxis_AFM[0] = std::stod(varValues[0]);
						mat.uniaxis_AFM[1] = std::stod(varValues[1]);
						mat.uniaxis_AFM[2] = std::stod(varValues[2]);
					}
					else if (varName == "AFM_EXCHANGE_ANISO") {
						mat.Aex_AFM1 = std::stod(varValues[0]);
						mat.Aex_AFM2 = std::stod(varValues[1]);
					}
					else if (varName == "AFM_EXCHANGE_COUPLING") mat.J_AFM = std::stod(varValues[0]);
					else if (varName == "AFM_LAMBDA_100") {
						mat.lamda100_AFM1 = std::stod(varValues[0]);
						mat.lamda100_AFM2 = std::stod(varValues[1]);
					}
					else if (varName == "AFM_LAMBDA_111") {
						mat.lamda111_AFM1 = std::stod(varValues[0]);
						mat.lamda111_AFM2 = std::stod(varValues[1]);
					}

					else if (varName == "IF_SPIN_PUMPING") {
						if (varValues[0] == "TRUE") mat.if_spin_pump = true;
						else mat.if_spin_pump = false;
					}
					else if (varName == "SPIN_HALL_ANGLE") mat.spin_hall_angle = std::stod(varValues[0]);
					else if (varName == "SPIN_MIXING_CONDUCTANCE") mat.spin_mix_cond = std::stod(varValues[0]);
					else if (varName == "SPIN_DIFFUSION_LENGTH") mat.spin_diffus_length = std::stod(varValues[0]);

					else if (varName == "IS_FERROELECTRIC") {
						if (varValues[0] == "TRUE") mat.if_FE = true;
						else mat.if_FE = false;
					}
					else if (varName == "FE_MASS") mat.FE_mass = std::stod(varValues[0]);
					else if (varName == "FE_DAMPING") mat.FE_damping = std::stod(varValues[0]);
					else if (varName == "FE_a1") mat.a1 = std::stod(varValues[0]);
					else if (varName == "FE_a11") mat.a11 = std::stod(varValues[0]);
					else if (varName == "FE_a12") mat.a12 = std::stod(varValues[0]);
					else if (varName == "FE_a111") mat.a111 = std::stod(varValues[0]);
					else if (varName == "FE_a112") mat.a112 = std::stod(varValues[0]);
					else if (varName == "FE_a123") mat.a123 = std::stod(varValues[0]);
					else if (varName == "FE_a1111") mat.a1111 = std::stod(varValues[0]);
					else if (varName == "FE_a1112") mat.a1112 = std::stod(varValues[0]);
					else if (varName == "FE_a1122") mat.a1122 = std::stod(varValues[0]);
					else if (varName == "FE_a1123") mat.a1123 = std::stod(varValues[0]);
					else if (varName == "FE_G11") mat.G11 = std::stod(varValues[0]);
					else if (varName == "FE_Q11") mat.Q11 = std::stod(varValues[0]);
					else if (varName == "FE_Q12") mat.Q12 = std::stod(varValues[0]);
					else if (varName == "FE_Q44") mat.Q44 = std::stod(varValues[0]);
					else if (varName == "FE_f11") mat.f11 = std::stod(varValues[0]);
					else if (varName == "FE_f12") mat.f12 = std::stod(varValues[0]);
					else if (varName == "FE_f44") mat.f44 = std::stod(varValues[0]);

					else if (varName == "RELATIVE_PERMITTIVITY") {
						if (varValues.size() != 9) {
							std::cerr << "Error in parsing relative permittivity (needs 9 components). Please check line number " << lineNumber << " in " << ssFileName << std::endl;
							return;
						}
						else {
							mat.r_permittivity11 = std::stod(varValues[0]);
							mat.r_permittivity12 = std::stod(varValues[1]);
							mat.r_permittivity13 = std::stod(varValues[2]);
							mat.r_permittivity21 = std::stod(varValues[3]);
							mat.r_permittivity22 = std::stod(varValues[4]);
							mat.r_permittivity23 = std::stod(varValues[5]);
							mat.r_permittivity31 = std::stod(varValues[6]);
							mat.r_permittivity32 = std::stod(varValues[7]);
							mat.r_permittivity33 = std::stod(varValues[8]);
						}
					}

					else if (varName == "ELECTRICAL_CONDUCTIVITY") {
						if (varValues.size() != 9) {
							std::cerr << "Error in parsing electrical conductivity (needs 9 components). Please check line number " << lineNumber << " in " << ssFileName << std::endl;
							return;
						}
						else {
							mat.conductivity11 = std::stod(varValues[0]);
							mat.conductivity12 = std::stod(varValues[1]);
							mat.conductivity13 = std::stod(varValues[2]);
							mat.conductivity21 = std::stod(varValues[3]);
							mat.conductivity22 = std::stod(varValues[4]);
							mat.conductivity23 = std::stod(varValues[5]);
							mat.conductivity31 = std::stod(varValues[6]);
							mat.conductivity32 = std::stod(varValues[7]);
							mat.conductivity33 = std::stod(varValues[8]);
						}
					}

					else if (varName == "PLASMA_OSCILLATION") {
						if (varValues.size() != 9) {
							std::cerr << "Error in parsing plasma oscillation (needs 9 values). Please check line number " << lineNumber << " in " << ssFileName << std::endl;
							return;
						}
						else {
							mat.comp_n1 = std::stod(varValues[0]);
							mat.omega_plasma_n1 = std::stod(varValues[1]);
							mat.tao_e_n1 = std::stod(varValues[2]);
							mat.comp_n2 = std::stod(varValues[3]);
							mat.omega_plasma_n2 = std::stod(varValues[4]);
							mat.tao_e_n2 = std::stod(varValues[5]);
							mat.comp_n3 = std::stod(varValues[6]);
							mat.omega_plasma_n3 = std::stod(varValues[7]);
							mat.tao_e_n3 = std::stod(varValues[8]);
						}
					}
				}
			}
		}
		
		mat.B1 = -3. / 2. * mat.lamda100 * (mat.c11 - mat.c12);
		mat.B2 = -3. * mat.lamda111 * mat.c44;

		mat.B1_AFM1 = -3. / 2. * mat.lamda100_AFM1 * (mat.c11 - mat.c12);
		mat.B1_AFM2 = -3. / 2. * mat.lamda100_AFM2 * (mat.c11 - mat.c12);
		mat.B2_AFM1 = -3. * mat.lamda111_AFM1 * mat.c44;
		mat.B2_AFM2 = -3. * mat.lamda111_AFM2 * mat.c44;

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

		check_material(mat);
		get_cijkl(mat);

		material_parameters[i] = mat;
	}
}

/**
 * Reads the geometry parameters from the "system_setting.in" file and sets the corresponding member variables.
 *
 * Potential updates:
 *	1. Rework input reading to use json or xml
 * 	2. User-specified input filename
 */
void geometry_parameters::readgeo() {

	const std::string ssFileName = "system_setting.in";

	uint32_t lineNumber = 0;

	std::ifstream ssFile(ssFileName);
	if (!ssFile.is_open()) {
		std::cerr << ssFileName << " could not be opened." << std::endl;
		std::cerr << "Check if the file exists." << std::endl;
		return;
	}

	while (!ssFile.eof()) {
		std::string line;
		std::string varName;

		std::getline(ssFile, line);
		lineNumber++;

		if (line.length() > 0) {
			if (line[0] == '#') {
				continue;
			}

			cleanUpLine(line, varName);

			std::vector<std::string> val;
			std::vector<std::string> varValues;

			std::stringstream ss(line);
			std::string token;
			while (std::getline(ss, token, ',')) {
				varValues.push_back(token);
			}

			if (varValues.size() == 0 || varValues[0] == "") {
				std::cerr << "Invalid input. Check parameter " << varName << " on line " << lineNumber << std::endl;
			}
			else {
				if (varName == "NX") {
					nx_system = std::stoi(varValues[0]);
					if (nx_system < 1)
						nx_system = 1;
				}
				else if (varName == "NY") {
					ny_system = std::stoi(varValues[0]);
					if (ny_system < 1)
						ny_system = 1;
				}
				else if (varName == "NZ") {
					nz_system = std::stoi(varValues[0]);
					if (nz_system < 1)
						nz_system = 1;
				}
				else if (varName == "DX") {
					dx = std::stod(varValues[0]);
				}
				else if (varName == "DY") {
					dy = std::stod(varValues[0]);
				}
				else if (varName == "DZ") {
					dz = std::stod(varValues[0]);
				}
				else if (varName == "PERIODIC") {
					if (varValues[0] == "TRUE") periodicX = true;
					else if (varValues[0] == "FALSE") periodicX = false;

					if (varValues[1] == "TRUE") periodicY = true;
					else if (varValues[1] == "FALSE") periodicY = false;

					if (varValues[2] == "TRUE") periodicZ = true;
					else if (varValues[2] == "FALSE") periodicZ = false;
				}
				else if (varName == "READ_STRUCT") {
					if (varValues[0] == "TRUE") if_read_struct = true;
					else if (varValues[0] == "FALSE") if_read_struct = false;
				}
				else if (varName == "NX_WORK") {
					nx_work = std::stoi(varValues[0]);
				}
				else if (varName == "NY_WORK") {
					ny_work = std::stoi(varValues[0]);
				}
				else if (varName == "SYSTEM_LAYERS") {
					num_layer_unit = std::stoi(varValues[0]);
					num_periods = std::stoi(varValues[1]);
					id_firstlayer = std::stoi(varValues[2]);
					id_lastlayer_unit = std::stoi(varValues[3]);
				}
				else if (varName == "NZ_LAYERS") {
					pt_nz_layer_unit = new unsigned int[num_layer_unit];

					for (long int i = 0; i < num_layer_unit; i++) {
						pt_nz_layer_unit[i] = std::stoi(varValues[i]);
					}
				}
				else if (varName == "PML") {
					if (varValues[0] == "TRUE") if_PML = true;
					else if (varValues[0] == "FALSE") if_PML = false;
				}
				else if (varName == "PML_DIRECTIONS") {
					if (true == if_PML) {

						if (varValues[0] == "TRUE") if_PML_Xs = true;
						else if (varValues[0] == "FALSE") if_PML_Xs = false;

						if (varValues[1] == "TRUE") if_PML_Xe = true;
						else if (varValues[1] == "FALSE") if_PML_Xe = false;

						if (varValues[2] == "TRUE") if_PML_Ys = true;
						else if (varValues[2] == "FALSE") if_PML_Ys = false;

						if (varValues[3] == "TRUE") if_PML_Ye = true;
						else if (varValues[3] == "FALSE") if_PML_Ye = false;

						if (varValues[4] == "TRUE") if_PML_Zs = true;
						else if (varValues[4] == "FALSE") if_PML_Zs = false;

						if (varValues[5] == "TRUE") if_PML_Ze = true;
						else if (varValues[5] == "FALSE") if_PML_Ze = false;
					}
					else {
						if_PML_Xs = false;
						if_PML_Xe = false;
						if_PML_Ys = false;
						if_PML_Ye = false;
						if_PML_Zs = false;
						if_PML_Ze = false;
					}
				}
				
				else if (varName == "PML_SIZE") {
					if (true == if_PML)
						PML_size = std::stoi(varValues[0]);
					else
						PML_size = 0;
				}
				else if (varName == "PML_MATERIAL_TYPE") {
					if (true == if_PML)
						PML_materialType = std::stoi(varValues[0]);
					else
						PML_materialType = 0;
				}
				else if (varName == "PML_KAPPA_MAX") {
					if (true == if_PML)
						kappaMax = std::stod(varValues[0]);
					else
						kappaMax = 0;
				}
				else if (varName == "PML_M") {
					if (true == if_PML)
						PML_m = std::stod(varValues[0]);
					else
						PML_m = 0;
				}
			}
		}
	}

	ssFile.close();
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
			material_cell(pt_geo->xS + x - 1, pt_geo->yS + y - 1, pt_geo->zS + z - 1) = id;
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

	logFile << "sigmaMax: " << sigmaMax << std::endl;
	logFile << "kappaMax: " << kappaMax << std::endl;
	logFile << "PML_d: " << PML_d << std::endl;
	logFile << "PML_m: " << PML_m << std::endl;
	logFile << "eta0: " << eta0 << std::endl;
	logFile << "maxReflErr: " << maxReflErr << std::endl;

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
	logFile << "PML_materialType: " << PML_materialType << std::endl;
	logFile << "PML boundaries - xS: " << xS << ", xE: " << xE << ", yS: " << yS << ", yE: " << yE << ", zS: " << zS << ", zE: " << zE << std::endl;

	// If there are other properties or arrays that need to be logged, you can add them here in a similar fashion.

	logFile.close();
}

