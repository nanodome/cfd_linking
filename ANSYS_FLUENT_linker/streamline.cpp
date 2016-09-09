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
