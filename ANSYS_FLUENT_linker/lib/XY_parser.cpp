/*
Original code by H2020 - NanoDome European Project (www.nanodome.eu, emanuele.ghedini@unibo.it)

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any
damages arising from the use of this software.

Permission is granted to anyone to use this software for any
purpose, including commercial applications, and to alter it and
redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software. If you use this
software in a product, an acknowledgment in the product documentation
would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and
must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source
distribution.
*/

#include "XY_parser.h"

// Auxiliary Functions
// Parse the particle list keeping trace of the time and the values
std::vector< std::vector<double> > getTime_Val(std::ifstream& file, std::vector< std::vector<double> >& time);

// Parse the particle list keeping trace of the values
std::vector< std::vector<double> > getVal(std::ifstream& file);

// Parse the file with the species
std::list<Species> getSpecies(std::ifstream& file);

// String Tokenizer
std::vector<std::string> tokenizerPar(const std::string& p_pcstStr, char delim);



XY_parser::XY_parser(){

	map_init();

}

void XY_parser::map_init(){

	// Init the map between string and enumeration for managing the properties file management
	prop_map["\"Static Pressure\""]				= Press;
	prop_map["\"Static Temperature\""]			= Temp;
	prop_map["\"X-Coordinate\""]				= X;
	prop_map["\"Y-Coordinate\""]				= Y;
	prop_map["\"Z-Coordinate\""]				= Z;
	prop_map["\"Scalar-2\""]					= Molar_C;
	prop_map["\"Species\""]						= P_Species;

	// List of properties names in XY format 
	prop_list = {"\"Static Pressure\"", "\"Static Temperature\"", "\"X-Coordinate\"", "\"Y-Coordinate\"", 
				 "\"Z-Coordinate\"", "\"Scalar-2\"", "\"Species\""}; // Scalar-2 is the dummy name given by FLUENT to the molar concentration
}


void XY_parser::parse(std::string filename, Raw_data& parsed_data){

	std::cout << "Parsing file: "<< filename<<std::endl;

	std::ifstream fluent_file(filename);
	std::string line;
	
	if (fluent_file.is_open()){
		// Get file header
		//first line;
		std::getline(fluent_file, line);

		// get property
		std::getline(fluent_file, line);

		// Identify the property
		std::string prop;
		size_t s_start;
		for (auto n_prop : prop_list){
			s_start = line.find(n_prop);
			if (s_start != -1){
				prop = line.substr(s_start, n_prop.size());
				break;
			}
		}

		// Skip space
		std::getline(fluent_file, line);
		
		if (s_start == -1)
		{
			std::cout << "Property not recognized, Quitting!!"<<std::endl;
			system("PAUSE");
			exit(0);
		}
		
		switch (prop_map[prop])
		{
		case Temp:
			std::cout << "Temperature File" << std::endl;
			parsed_data.temperatures = getTime_Val(fluent_file, parsed_data.times);
			break;

		case Press:
			std::cout << "Pressure File" << std::endl;
			parsed_data.pressures = getVal(fluent_file);
			break;

		case X:
			std::cout << "X File" << std::endl;
			parsed_data.x = getVal(fluent_file);
			break;

		case Y:
			std::cout << "Y File" << std::endl;
			parsed_data.y = getVal(fluent_file);
			break;

		case Z:
			std::cout << "Z File" << std::endl;
			parsed_data.z = getVal(fluent_file);
			break;

		case Molar_C:
			std::cout << "Molar Concentration File" << std::endl;
			parsed_data.molar_cs = getVal(fluent_file);
			break;

		case P_Species:
			std::cout << "Species File" << std::endl;
			parsed_data.species = getSpecies(fluent_file);
			break;

		default: // if property not recognized -> QUITTING
			std::cout << "Unrecognized  File" << std::endl;
			system("PAUSE");
			exit(0);
			break;
		}

	}
	else{ // if file is not opened
		std::cout << "file: " << filename << " is not open" << std::endl;
	}

	fluent_file.close();

}


// Parse the particle list keeping trace of the time and the values
std::vector< std::vector<double> > getTime_Val(std::ifstream& file, std::vector< std::vector<double> >& time){

	std::vector< std::vector<double> > res;

	std::string line;
	std::vector<double> tmp;
	std::vector<double> t_tmp;


	while (std::getline(file, line)){
		if (line.compare("") == 0){
			
			continue;
		}
		else if (line.find("particle") != std::string::npos){
			
			continue;
		}
		else if (line.compare(")") == 0){
			
			res.push_back(tmp);
			time.push_back(t_tmp);

			tmp.clear();
			t_tmp.clear();

		}
		else{
			std::vector<std::string> tokens = tokenizerPar(line, '\t');
			double time_sample = atof(tokens[0].c_str());
			double val = atof(tokens[1].c_str());
			t_tmp.push_back(time_sample);
			tmp.push_back(val);
		}
	}

	return res;
}

// Parse the particle list keeping trace of the values
std::vector< std::vector<double> > getVal(std::ifstream& file){

	std::vector< std::vector<double> > res;

	std::string line;
	std::vector<double> tmp;
	
	

	while (std::getline(file, line)){
		if (line.compare("") == 0){
			continue;
		}
		else if (line.find("particle") != std::string::npos){
			continue;
		}
		else if (line.compare(")") == 0){
			
			res.push_back(tmp);
			tmp.clear();
			continue;
		}
		else{
			std::vector<std::string> tokens = tokenizerPar(line, '\t');
			double val = atof(tokens[1].c_str());
			tmp.push_back(val);
		}
	}

	return res;
}

std::list<Species> getSpecies(std::ifstream& file){

	std::list<Species> res;

	std::string line;

	// Skip Space
	std::getline(file, line);

	// Get species line
	std::vector<std::string> tokens = tokenizerPar(line, ' ');

	for (int s = 0; s < tokens.size(); s++){
		Species tmp(tokens[s]);

		res.push_back(tmp);
	}

	return res;


}

std::vector<std::string> tokenizerPar(const std::string& p_pcstStr, char delim)
{
    std::vector<std::string> tokens;
    std::stringstream mySstream(p_pcstStr);
    std::string temp;

	while (getline(mySstream, temp, delim)){
		if (!temp.empty())
		tokens.push_back(temp);
	}

    return tokens;
}




