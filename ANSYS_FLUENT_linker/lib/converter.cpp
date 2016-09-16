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

#include "converter.h"

Converter::Converter(std::string _dir_path, std::string _ret_path){

	dir_path = _dir_path;

	ret_path = _ret_path;
}

void Converter::Convert(std::string choice){

	// get filenames from the directory
	get_files();

	// declare parser for analyzing XY fluent files
	XY_parser File_parser;

	// Parse all files in the directory
	for (auto file_n : xy_files){
		std::string path = dir_path + "/" + file_n;
		File_parser.parse(path, data);
	}
	
	// Create Streamlines
	for (int s = 0; s < data.get_n_streams(); s++){

		streamline s_tmp(s, data.get_temperatures(s), data.get_pressures(s), data.get_times(s), 
						 data.get_species(), data.get_molar_cs(s), 
						 data.get_x(s), data.get_y(s), data.get_z(s));

		streams.push_back(s_tmp);
	}

	// Write XML 

	if (choice.compare("x") == 0){ // Print only XML

		streamlinefileXML XML_file(ret_path + "/NanoDomeStreamlines.xml", streams);
		XML_file.write_Streamlines(ret_path + "/NanoDomeStreamlines.xml");
	}
	else if (choice.compare("j") == 0){ // Pront only JSON

		streamlinefileJSON JSON_out(ret_path + "/NanoDomeStreamlines.json");
		JSON_out.write_streamlines(streams);
	}
	else{ // Print XML and JSON

		streamlinefileXML XML_file(ret_path + "/NanoDomeStreamlines.xml", streams);
		XML_file.write_Streamlines(ret_path + "/NanoDomeStreamlines.xml");

		streamlinefileJSON JSON_out(ret_path + "/NanoDomeStreamlines.json");
		JSON_out.write_streamlines(streams);
	}
	
	return;

}

void Converter::get_files(){

	DIR *folder = opendir(dir_path.c_str());
	
	std::string filename;

	dirent *content = readdir(folder);
	
	while (content != NULL){
		filename = content->d_name;
		if (filename.compare(".") != 0 && filename.compare("..") != 0){
			xy_files.push_back(filename);
			
		}
		content = readdir(folder);
	}

	closedir(folder);

}