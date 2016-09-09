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