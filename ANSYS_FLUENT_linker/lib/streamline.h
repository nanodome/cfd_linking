#ifndef STREAMLINE
#define STREAMLINE

#include <list>
#include <string>
#include <vector>
#include <iostream>

#include "species.h"

class streamline{

    /// Streamline ID
    int ID;

    /// Temperatures of the streamline in function of time [K]
    std::vector<double> T;

    /// Pressures of the streamline in funcion of time [Pa]
    std::vector<double> P;

    /// Time samples of the streamline [sec.]
    std::vector<double> Time;

    /// Species in the streamline in funtion of time
    std::list<Species> species;

    /// Species molar concentration for each time sample
    std::vector<double> Molar_Conc;

	std::vector<double> X;

	std::vector<double> Y;

	std::vector<double> Z;

public:

    /// Default Constructor
    streamline();

    /// Parametric Constructor
    streamline(int _ID, std::vector<double> _T, std::vector<double> _P,
               std::vector<double> _Time, std::list<Species> _species,
			   std::vector<double> _Molar_Conc, std::vector<double> _X, std::vector<double> _Y, std::vector<double> _Z);

    /// Return Temperature [K] array
    std::vector<double> get_Temp() const {return T;}

    /// Return Pressure [Pa] array
    std::vector<double> get_Press() const {return P;}

    /// Return Times [sec]
    std::vector<double>& get_Time() {return Time;}

	/// Return Molar Concentration for each species [sec]
	std::vector<double> get_Molar() const { return Molar_Conc; }

    /// Return Species
    std::list<Species> get_Species() const {return species;}

	/// Return Streamline ID
	int get_ID() const { return ID; }

	/// Return X Positions
	std::vector<double> get_x() const { return X; }

	/// Return Y Positions
	std::vector<double> get_y() const { return Y; }

	/// Return Z Positions
	std::vector<double> get_z() const { return Z; }

    /// Print Streamline
    void printstream();




};

#endif // STREAMLINE

