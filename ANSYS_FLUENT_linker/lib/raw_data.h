#ifndef RAW_DATA_H
#define RAW_DATA_H

#include <vector>
#include <list>
#include <string>

#include "species.h"

class Raw_data{

public:

	/// raw data from the times extracted from XY file from ANSYS FLUENT
	std::vector< std::vector< double > > times;

	/// raw data from the tenperatures extracted from XY file from ANSYS FLUENT
	std::vector< std::vector< double > > temperatures;

	/// raw data from the pressures extracted from XY file from ANSYS FLUENT
	std::vector< std::vector< double > > pressures;

	/// raw data from the species molar concentrations extracted from XY file from ANSYS FLUENT
	std::vector< std::vector< double > > molar_cs;

	/// raw data from the x coordinates extracted from XY file from ANSYS FLUENT
	std::vector< std::vector< double > > x;

	/// raw data from the y coordinates extracted from XY file from ANSYS FLUENT
	std::vector< std::vector< double > > y;

	/// raw data from the z coordinates extracted from XY file from ANSYS FLUENT
	std::vector< std::vector< double > > z;

	/// raw data from the species extracted from XY file from ANSYS FLUENT
	std::list< Species > species;


	/// Blank Contructor
	Raw_data();

	/*void update_times(std::vector< double > time_v)			{ times.push_back(time_v); }

	void update_temperatures(std::vector< double > temp_v)	{ temperatures.push_back(temp_v); }

	void update_pressures(std::vector< double > press_v)	{ pressures.push_back(press_v); }

	void update_x(std::vector< double > x_v)				{ x.push_back(x_v); }
		
	void update_y(std::vector< double > y_v)				{ y.push_back(y_v); }

	void update_z(std::vector< double > z_v)				{ z.push_back(z_v); }

	void update_molar_c(std::vector< double > molar_v)		{ molar_cs.push_back(molar_v); }*/

	/// Returns the number of streamlines tracked
	int get_n_streams() const{return times.size();}

	/// Initiate the time samples array
	/// \param std::vector< std::vector< double > > time_v: time samples extracted from XY files
	void init_times(std::vector< std::vector< double > >& time_v)			{ times = time_v; }

	/// Initiate the temperatures samples arrays
	/// \param std::vector< std::vector< double > > temp_v: temperatures samples extracted from XY files
	void init_temperatures(std::vector< std::vector< double > >& temp_v)	{ temperatures = temp_v; }

	/// Initiate the pressures samples arrays
	/// \param std::vector< std::vector< double > >  press_v: pressures samples extracted from XY files
	void init_pressures(std::vector< std::vector< double > >& press_v)		{ pressures = press_v; }

	/// Initiate the x coordinates arrays
	/// \param std::vector< std::vector< double > >  x_v: x coordinates extracted from XY files
	void init_x(std::vector< std::vector< double > >& x_v)					{ x = x_v; }
	
	/// Initiate the y coordinates arrays
	/// \param std::vector< std::vector< double > >  y_v: y coordinates extracted from XY files
	void init_y(std::vector< std::vector< double > >& y_v)					{ y = y_v; }

	/// Initiate the z coordinates arrays
	/// \param std::vector< std::vector< double > >  z_v: z coordinates extracted from XY files
	void init_z(std::vector< std::vector< double > >& z_v)					{ z = z_v; }

	/// Initiate the molar concentration samples arrays
	/// \param std::vector< std::vector< double > >  molar_v: molar concentration for each species extracted from XY files
	void init_molar_c(std::vector< std::vector< double > >& molar_v)		{ molar_cs = molar_v; }



	std::vector<double>		get_times(int idx) const {std::vector< double > dummy = {0};
													  return (times.size() != 0) ? times[idx] : dummy; };

	std::vector< double >	get_temperatures(int idx) const {std::vector< double > dummy = {0};
															 return (temperatures.size() != 0) ? temperatures[idx] : dummy;};

	std::vector< double >	get_molar_cs(int idx) const {std::vector< double > dummy = {0}; 
														 return (molar_cs.size() != 0) ? molar_cs[idx] : dummy;};

	std::vector< double >	get_x(int idx) const {std::vector< double > dummy = {0}; 
												  return (x.size() != 0) ? x[idx] : dummy;};

	std::vector< double >	get_y(int idx) const {std::vector< double > dummy = {0}; 
												  return (y.size() != 0) ? y[idx] : dummy;};

	std::vector< double >	get_z(int idx) const {std::vector< double > dummy = {0}; 
												  return (z.size() != 0) ? z[idx] : dummy;};

	std::vector< double >	get_pressures(int idx) const {std::vector< double > dummy = {0}; 
														  return (pressures.size() != 0) ? pressures[idx] : dummy;};

	std::list< Species >	get_species() const {return species;};




};

#endif