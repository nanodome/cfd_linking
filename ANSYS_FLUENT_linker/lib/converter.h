/*
Original code by H2020 - NanoDome European Project (www.nanodome.eu, https://github.com/nanodome/cfd_linking.git)

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

#ifndef CONVERTER_H
#define CONVERTER_H


#include <vector>
#include <string>

// Preprocessor defines for using the libraries relative to the folder management 
#if (defined(_WIN32) || defined(WIN32))  && !defined(__GNUC__)
#include "win_dir_libs/dirent.h"
#else
#include <dirent.h>
#endif

#include "streamline.h"
#include "XY_parser.h"
#include "raw_data.h"
#include "streamlinefileXML.h"
#include "streamlinefileJSON.h"

class Converter{

	/// Path of the directory containing ANSYS FLUENT .xy files
	std::string dir_path;

	/// Path of the directories containing the returned file
	std::string ret_path;

	/// List of the names of the files
	std::list<std::string> xy_files;

	/// streamlines extractred from the files
	std::vector<streamline> streams;

	/// Raw data parsed from the files
	Raw_data data;

	/// performs 
	void get_files();

public:

	/// Contructor
	/// \param std::string _dir_path: directory path
	Converter(std::string _dir_path, std::string _ret_path);

	/// Converts XY ANSYS FLUENT files to XML or JSON NanoDome format
	void Convert(std::string choice);

};


#endif