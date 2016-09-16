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