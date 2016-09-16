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

#ifndef STREAMLINEFILEJSON_H
#define STREAMLINEFILEJSON_H

/*IMPORTANT!!! when you change the platform and the compiler remember that the preprocessor definition in JSON.h 
               can give problems */


#include "JSONFile.h"
#include "streamline.h"



class streamlinefileJSON :public JSONfile{

	

public:

	/// Contructor
	/// \param std::string _file: path to the JSON file
	streamlinefileJSON(std::string _file);

	/// Parse the JSON file
	std::vector<streamline> parse();

	/// Write the JSON file starting from a set of streamlines
	void write_streamlines(std::vector<streamline>& _streams);

};



#endif // STREAMLINEFILEJSON_H