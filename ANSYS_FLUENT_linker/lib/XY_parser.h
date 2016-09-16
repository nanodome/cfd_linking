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