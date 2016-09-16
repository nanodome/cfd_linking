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