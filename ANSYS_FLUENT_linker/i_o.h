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

