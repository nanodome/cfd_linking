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

#include <fstream>
#include <sstream>
#include <cstring>
#include <list>

#include "species.h"
#include "streamlinefileJSON.h"



streamlinefileJSON::streamlinefileJSON(std::string _file){

	JSON_file_path = _file;
}

std::vector<streamline> streamlinefileJSON::parse(){

	std::vector<streamline> streams;

	std::string line;
	std::ifstream j_file(JSON_file_path.c_str());
	std::stringstream buffer;
	buffer << j_file.rdbuf();
	line = buffer.str();

	int tmp_n_streams;
	double tmp_start_t, tmp_end_t;

	j_file.close();

	JSONValue* j_val = JSON::Parse(line.c_str());
	
	root = j_val->AsObject();
	
	read_JSON_value(root, "N_stream" , tmp_n_streams);

	// Get start time and end time
	read_JSON_value(root, "T_star", tmp_start_t);
	read_JSON_value(root, "T_end", tmp_end_t);

	// Get the streamlines data
	if (root.find(L"Streams") != root.end() && root[L"Streams"]->IsArray()){

		// Get each stremline data
		JSONArray s_array = root[L"Streams"]->AsArray();

		for (int s = 0; s < s_array.size(); s++){
			JSONObject stream = s_array[s]->AsObject();

			int s_id; int n_t_samp, n_species;
			std::vector<double> s_Time, s_T, s_P, s_X, s_Y, s_Z, s_Con;
			std::list<Species> s_Species;

			// get stremline index
			read_JSON_value(stream, "Id", s_id);

			//get streamline number of time samples
			read_JSON_value(stream, "N_time_samples", n_t_samp);

			//get timesamples
			read_JSON_array(stream, "Time_samples", s_Time);

			// get Temperaturees
			read_JSON_array(stream, "T", s_T);

			// get pressures
			read_JSON_array(stream, "P", s_P);

			// get X
			read_JSON_array(stream, "X", s_X);

			// get Y
			read_JSON_array(stream, "Y", s_Y);

			// get Z
			read_JSON_array(stream, "Z", s_Z);

			// get number of species
			read_JSON_value(stream, "N_species", n_species);

			// Get Species
			if (stream.find(L"Species") != stream.end() && stream[L"Species"]->IsArray()){

				JSONArray arr = stream[L"Species"]->AsArray();
				for (int i = 0; i < arr.size(); i++){
					std::wstring formula = arr[i]->AsString();
					std::string s_species(formula.begin(), formula.end());
					Species s(s_species);
					s_Species.push_back(s);
				}
			}

			// Get molar
			read_JSON_array(stream, "Molar_c", s_Con);

			// Create the streamline
			streamline j_stream(s_id, s_T, s_P, s_Time, s_Species, s_Con, s_X, s_Y, s_Z);

			// Add the streamline to the set
			streams.push_back(j_stream);
		}
	}

	// Free root
	root.clear();

	return streams;

}

void streamlinefileJSON::write_streamlines(std::vector<streamline>& _streams){

	std::cout << "Creating JSON File: "<< JSON_file_path <<std::endl;

	// Output File
	std::wofstream out_file((JSON_file_path).c_str(), std::wfstream::out);
	
	// Number of streams
	int n_streams = _streams.size();
	
	root[L"N_stream"] = new JSONValue((int)n_streams);

	
	// Starting time
	root[L"T_start"] = new JSONValue((double)_streams[0].get_Time()[0]);

	// Ending Time
	root[L"T_end"] = new JSONValue((double)_streams[0].get_Time()[_streams[0].get_Time().size() - 1]);

	// Array of streamlines JSON Objects
	JSONArray j_streams_arr;

	// Insert each streamline
	for (auto s_l: _streams){

		// Create streamline JSON object
		JSONObject j_stream;

		// Add stream id
		j_stream[L"Id"] = new JSONValue((int)s_l.get_ID());

		// Add stream number of time samples
		j_stream[L"N_time_samples"] = new JSONValue((int)s_l.get_Time().size());
		
		// Add Time Samples
		std::vector<double> v_Time = s_l.get_Time();
		j_stream[L"Time_samples"] = new JSONValue(get_JSON_array(v_Time, v_Time.size()));

		// Add Temperatures
		std::vector<double> v_T  = s_l.get_Temp();
		j_stream[L"T"] = new JSONValue(get_JSON_array(v_T, v_T.size()));

		// Add Pressures
		std::vector<double> v_P = s_l.get_Press();
		j_stream[L"P"] = new JSONValue(get_JSON_array(v_P, v_P.size()));


		// Retrieve species
		std::list<Species> s_list = s_l.get_Species();
		
		// Add Species Cardinality
		j_stream[L"N_species"] = new JSONValue((int)s_list.size());

		// Add Species formula
		JSONArray species_arr;
		for (auto sp: s_list){
			std::string formula = sp.get_formula();
			std::wstring w_species(formula.begin(), formula.end());
			species_arr.push_back(new JSONValue(w_species));
		}
		j_stream[L"Species"] = new JSONValue(species_arr);

		// Add Molar Concentration 

		std::vector<double> v_Mc = s_l.get_Molar();
		j_stream[L"Molar_c"] = new JSONValue(get_JSON_array(v_Mc, v_Mc.size()));

		// Add X
		std::vector<double> v_X = s_l.get_x();
		j_stream[L"X"] = new JSONValue(get_JSON_array(v_X, v_X.size()));

		// Add Y
		std::vector<double> v_Y = s_l.get_y();
		j_stream[L"Y"] = new JSONValue(get_JSON_array(v_Y, v_Y.size()));

		// Add Z
		std::vector<double> v_Z = s_l.get_z();
		j_stream[L"Z"] = new JSONValue(get_JSON_array(v_Z, v_Z.size()));

		// Insert the streamline in the array
		j_streams_arr.push_back(new JSONValue(j_stream));

	}

	// Insert the streamlines inside the main object
	root[L"Streams"] = new JSONValue(j_streams_arr);
	
	// Create the final object
	JSONValue* out_value = new JSONValue(root);

	// Create the final string
	std::wstring out_string = out_value->Stringify();

	// Print it on file
	out_file << out_string;

	// Close the file
	out_file.close();

	// Free memory
	delete out_value;

	return;

}


