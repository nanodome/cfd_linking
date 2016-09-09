#ifndef SPECIES_H
#define SPECIES_H

#include "nanodome.h"

#include <string>
#include <cmath>


class Species {

// fundamental properties
    std::string formula; ///< species formula (e.g. "ZnO")
    std::string name; ///< iupac species name (e.g. "zinc monoxide")

    double mass; ///< species mass [kg]

    double T_melt; ///< melting temperature [K]

    double sigma; ///< L-J sigma value [m]
    double eps; ///< L-J epsilon value [J]

// CNT oriented properties

    // densities
    double bulk_density_liq; ///< density of the liquid bulk phase [kg/m3]
    double bulk_density_sol; ///< density of the solid bulk phase [kg/m3]

    // surface tension is expressed as sten_A - sten_B*(T - sten_C) [N/m]
    double s_ten_A; ///< surface tension coefficient A [N/m]
    double s_ten_B; ///< surface tension coefficient B [N/(m K)]
    double s_ten_C; ///< surface tension coefficient C [K]

    // saturation pressure is expressed as log10(psat) = (psat_A - psat_B/T) * 1.01e5 [Pa]
    double p_sat_A; ///< saturation pressure coefficient A
    double p_sat_B; ///< saturation pressure coefficient B [K]

public:

    Species(std::string _formula);

    ///< get formula
    std::string get_formula() const { return formula; }

    ///< species mass [kg]
    double get_mass() const { return mass; }

    ///< species L-J sigma [m]
    double get_sigma() const { return sigma; }

    ///< species L-J epsilon [J]
    double get_epsilon() const { return eps; }

    ///< bulk material surface tension [N/m]
    /// \param T temperature [K]
    double s_ten(double T) const { return s_ten_A-s_ten_B*(T-s_ten_C); }

    ///< saturation pressure [Pa]
    /// \param T temperature [K]
    double p_sat(double T) const { return 1.01e5 * pow(10.0,(p_sat_A-(p_sat_B/T))); }

    ///< saturation density [#/m3]
    /// \param T temperature [K]
    double n_sat(double T) const { return p_sat(T)/(K_BOL*T); }

    ///< molecular volume based on liquid density [m3]
    double m_volume() const { return mass/(AMU*N_AVO*1000*bulk_density_liq); }

    ///< molecular surface based on liquid density and spherical assumption [m2]
    double m_surface() const { return 4.8360*pow(m_volume(),2./3.); }

    double get_bulk_density(double T) const { return (T<T_melt) ? bulk_density_sol : bulk_density_liq; }
};

#endif // SPECIES_H
