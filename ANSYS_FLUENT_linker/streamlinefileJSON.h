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