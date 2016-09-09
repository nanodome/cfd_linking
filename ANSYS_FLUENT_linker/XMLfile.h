#ifndef _XMLFILE_H
#define _XMLFILE_H

#include "i_o.h"
#include "tinyxml2/tinyxml2.h"

using namespace tinyxml2;

class XMLfile :public i_o{

protected:
	/// XML file Path
	std::string path;

	/// XML file (from tinyxml)
	XMLDocument doc;
	
	/// Create Node for the XML File
	/// \param XMLElement *cursor_ptr: pointer for expanding the <TAG> sub_tree
	/// \param std::string TAG: Name of the node
	/// \param string: Value of the node
	void Create_XML_Node(XMLElement *root_e, std::string TAG, std::string text);

	/// Check if the requested XML TAG is correct
	/// \param XMLElement *e_ptr: Pointer to the tinyxml ElementXML to check
	/// \param std::string TAG
	bool check_Tag(XMLElement *e_ptr, std::string TAG);

	/// Print out a detalled error and stops the execution
	/// bool status: status varible from other processes
	/// std::string TAG: XML tag giving problems
	void Error_Check(bool status, std::string TAG);

};

#endif