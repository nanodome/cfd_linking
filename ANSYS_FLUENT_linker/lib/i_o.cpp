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

#include "i_o.h"

i_o::i_o(){}


void i_o::log_entry(std::string entry){



}


void i_o::error_entry_not_blocking(std::string _error){

	std::cout << "ERROR!! Check errors.dat\n";
	freopen(err_file_path.c_str(), "a", stderr);
	std::cerr << _error << std::endl;
	fclose(stderr);
}

void i_o::error_entry_blocking(std::string _error){

	std::cout << "ERROR!! Check errors.dat\n";
	freopen(err_file_path.c_str(), "a", stderr);
	std::cerr << _error << std::endl;
	fclose(stderr);
	system("PAUSE");
	exit(0);

}

