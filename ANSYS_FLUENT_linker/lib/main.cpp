#include <string>

#include "converter.h"

int main(int argc, char *argv[]){
	
	std::string in_path, out_path, choice;

	/*------ FOR TESTING PURPOSES ----------*/
	//in_path = "fluent";
	//out_path = "test";
	//choice = "xj"; // x j xj

	in_path = argv[1];
	out_path = argv[2];
	choice = argv[3];

	Converter XYtoXML(in_path, out_path);
	
	XYtoXML.Convert(choice);

	//system("PAUSE");

	return 0;
}