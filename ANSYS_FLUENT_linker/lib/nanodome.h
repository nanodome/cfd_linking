/*
Original code by H2020 - NanoDome European Project (www.nanodome.eu, https://github.com/nanodome/cfd_linking.git)

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

#ifndef NANODOME
#define NANODOME


// useful constants
#define K_BOL 1.380650524e-23    // [J/K]
#define E_CHARGE 1.602176565e-19 // [C]
#define AMU 1.660538921e-27      // [kg]
#define N_AVO 6.0221408571e23    // [#/mol]


// normalization units
#define MASS_UNIT 1.66053892e-027   //kg (= 1 amu)
#define LENGTH_UNIT 1.0e-010        //m (= 1 angstrom)
#define TIME_UNIT 1.01805055e-014   //s
#define ENERGY_UNIT 1.60217657e-019 //J (= 1 ev)
#define FORCE_UNIT 1.60217657e-009  //N
#define TEMPERATURE_UNIT (1.0/8.617332478e-5) //K


// useful math macros
#define sqr(x) ((x)*(x))


// activate event log
#define VERBOSE

// parameters used for SHAKE algorithm
#define SHAKE_CONV_CRITERION 1e-4
#define SHAKE_MAX_ITERATIONS  100// default 1000

#endif // NANODOME

