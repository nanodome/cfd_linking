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

#include "XMLfile.h"

void XMLfile::Create_XML_Node(XMLElement *root_e, std::string TAG, std::string text){

	XMLElement *p_new_elem = doc.NewElement(TAG.c_str());
	p_new_elem->SetText(text.c_str());
	root_e->InsertEndChild(p_new_elem);

}

bool XMLfile::check_Tag(XMLElement *e_ptr, std::string TAG){

	bool check = false;

	if (e_ptr == nullptr || strcmp(e_ptr->Name(), TAG.c_str()))
	{
		std::string error = "Tag " + TAG + " NOT FOUND" + "\n";
		error_entry_not_blocking(error);
		return check;
	}
	else{
		check = true;
		std::cout << "Tag: " << TAG << " FOUND" << std::endl;
		return check;
	}

}

void XMLfile::Error_Check(bool status, std::string TAG){

	if (!status){

		std::string error = "file: " + path + " Streamlines XML file Parsing Error at tag: " + TAG + "\n";
		error_entry_blocking(error);
	}
}