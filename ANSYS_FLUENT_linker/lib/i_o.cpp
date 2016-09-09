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

