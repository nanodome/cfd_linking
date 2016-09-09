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