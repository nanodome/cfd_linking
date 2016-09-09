#ifndef XY_PARSER_H
#define XY_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "raw_data.h"

class XY_parser{

	
	/// Enumerator for files properties
    enum stream_prop {Temp, Press, X, Y, Z, Molar_C, P_Species};

    /// map for Strings and properties
    std::map<std::string, stream_prop> prop_map;

	/// vector of possible properies
	std::vector<std::string> prop_list;

public:

	/// Object blank constructor
	XY_parser();

	/// Parse the given XY file
	/// \param std::string filename: file name
	/// \param raw_data& parsed_data reference to the raw data extracted from the files
	void parse(std::string filename, Raw_data& parsed_data);

	/// Inititialize support data for the file parsing
	void map_init();

};



#endif