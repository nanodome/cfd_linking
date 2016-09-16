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

#include "streamline.h"

streamline::streamline(){}

streamline::streamline(int _ID, std::vector<double> _T, std::vector<double> _P,
           std::vector<double> _Time, std::list<Species> _species,
		   std::vector<double> _Molar_Conc, std::vector<double> _X, std::vector<double> _Y, std::vector<double> _Z){

    ID = _ID;
    T = _T;
    P = _P;
    Time = _Time;
    species = _species;
    Molar_Conc = _Molar_Conc;
	X = _X;
	Y = _Y;
	Z = _Z;
}

void streamline::printstream(){

    std::cout<<"ID: "<< ID <<std::endl;

    std::cout<<"Time Samples"<<std::endl;
    for(auto t = Time.begin(); t != Time.end(); ++t) std::cout<<(*t) << " ";
	std::cout << std::endl;

    std::cout<<"Temperatures"<<std::endl;
    for(auto t = T.begin(); t != T.end(); ++t) std::cout<<(*t) << " ";
	std::cout << std::endl;

    std::cout<<"Pressure"<<std::endl;
    for(auto t = P.begin(); t != P.end(); ++t) std::cout<<(*t) << " ";
	std::cout << std::endl;

    std::cout<<"Species"<<std::endl;
    for(auto t = species.begin(); t != species.end(); ++t)
        std::cout<<(*t).get_formula() << " ";
	std::cout << std::endl;

    std::cout<<"Molar Concentration"<<std::endl;
    for(auto t = Molar_Conc.begin(); t != Molar_Conc.end(); ++t) std::cout<<(*t) << " ";
	std::cout << std::endl;


    std::cout<<std::endl;


}
