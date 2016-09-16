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

#ifndef _JSONFILE_H
#define _JSONFILE_H


#include <vector>

#include "JSON/JSON.h"
#include "JSON/JSONValue.h"

#include "i_o.h"

class JSONfile :public i_o{

protected:
	
	/// path to the JSON file
	std::string JSON_file_path;

	// root of the JSON file
	JSONObject root;


	/// Template class for creating JSON array from std array (Implemented here because are template)
	/// \param T& vect: vector to convert (numerical)
	/// \param int dim: vector lenght
	template<class T>
	JSONArray get_JSON_array(std::vector<T>& vect, int dim){
		JSONArray arr;

		for (int i = 0; i < dim; i++){
			arr.push_back(new JSONValue(vect[i]));
		}

		return arr;
	}

	/// Template class for managing JSON values retrieving from file (Implemented here because are template)
	/// \param JSONObject j_obj: inspected JSON Object
	/// \param std::string tag: tag searched
	/// \param T& val: return value (numerical, int or float)
	template<class T>
	void read_JSON_value(JSONObject j_obj, std::string tag, T& val){
		std::wstring w_tag(tag.begin(), tag.end());

		if (j_obj.find(w_tag) != j_obj.end() && j_obj[w_tag]->IsNumber()){

			val = j_obj[w_tag]->AsNumber();
		}
		else{
			std::string error = "JSON Value -" + tag + "- Not found at file: " + JSON_file_path + "\n";
			error_entry_blocking(error);
		}
	}

	/// Template class for managing JSON arrays retrieving from file (Implemented here because are template)
	/// \param JSONObject j_obj: inspected JSON Object
	/// \param std::string tag: tag searched
	/// \param T& val: return vector (numerical, int or float)
	template<class T>
	void read_JSON_array(JSONObject j_obj, std::string tag, std::vector<T>& val){

		std::wstring w_tag(tag.begin(), tag.end());

		if (j_obj.find(w_tag) != j_obj.end() && j_obj[w_tag]->IsArray()){

			JSONArray arr = j_obj[w_tag]->AsArray();
			for (int i = 0; i < arr.size(); i++)
				val.push_back(arr[i]->AsNumber());

		}
		else{
			std::string error = "JSON array " + tag + "Not found at file: " + JSON_file_path + "\n";
			error_entry_blocking(error);
		}
	}

};

#endif