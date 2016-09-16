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

#ifndef I_O_H
#define I_O_H

#include <iostream>
#include <string>


class i_o{

	/// errors file path
	std::string err_file_path = "err/errors.dat";

	/// log file path
	std::string log_file_path = "log/log.dat";

public:

	/// Contructor
	i_o();

	/// function for writing to the log file
	/// \param std::string entry: entry to write
	void log_entry(std::string _entry);

	/// function for writing to the error output (blocking)
	/// \param std::string error: error description
	void error_entry_not_blocking(std::string _error);

	/// function for writing to the error output (blocking)
	/// \param std::string error: error description
	void error_entry_blocking(std::string _error);


};

#endif // I_O

