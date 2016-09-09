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