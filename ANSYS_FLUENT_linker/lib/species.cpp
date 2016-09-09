#include "species.h"

#include <iostream>


// for now only Si and Ar are implemented
Species::Species(std::string _formula) {

	if (_formula == "Si") {

		formula = "Si";
		name = "silicon";

		mass = 28.085*AMU;

		T_melt = 1687.;

		sigma = 3.30e-10;
		eps = 4.37e-20;

		bulk_density_liq = 2570.;
		bulk_density_sol = 2329.;

		s_ten_A = 0.732;
		s_ten_B = 0.000086;
		s_ten_C = 1685.;

		p_sat_A = 7.5341;
		p_sat_B = 23399.;

	}
	else if (_formula == "Ar") {

		formula = "Ar";
		name = "argon";

		mass = 39.948*AMU;

		T_melt = 83.81;

		sigma = 3.41e-10;
		eps = 1.65e-21;

		bulk_density_liq = 1395.4;
		bulk_density_sol = 1395.4;

		s_ten_A = 0.;
		s_ten_B = 0.;
		s_ten_C = 0.;

		p_sat_A = 0.;
		p_sat_B = 0.;

	}
	else {
		if (_formula == "Fe"){

			formula = "Fe";
			name = "Iron";

			mass = 39.948*AMU;

			T_melt = 83.81;

			sigma = 3.41e-10;
			eps = 1.65e-21;

			bulk_density_liq = 1395.4;
			bulk_density_sol = 1395.4;

			s_ten_A = 0.;
			s_ten_B = 0.;
			s_ten_C = 0.;

			p_sat_A = 0.;
			p_sat_B = 0.;
		}
		else{

			std::cout << "Species " << _formula << " does not exist!!!" << std::endl;

			exit(1);
		}
	}
}
