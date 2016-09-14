/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::streamParcelInjectionData

Description
    Container class to provide injection data for stream parcels

SourceFiles
    streamParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef streamParcelInjectionData_H
#define streamParcelInjectionData_H

#include "dictionary.H"
#include "vector.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class streamParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const streamParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    streamParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
               Class streamParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class streamParcelInjectionData
{
protected:

    // Parcel properties

        //- Position [m]
        point x_;

        //- Velocity [m/s]
        vector U_;

        //- Diameter [m]
        scalar d_;

        //- Density [kg/m3]
        scalar rho_;

        //- Mass flow rate [kg/s]
        scalar mDot_;


public:

    //- Runtime type information
    TypeName("streamParcelInjectionData");

    // Constructors

        //- Null constructor
        streamParcelInjectionData();

        //- Construct from dictionary
        streamParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        streamParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~streamParcelInjectionData();


    // Access

        //- Return const access to the position
        inline const point& x() const;

        //- Return const access to the velocity
        inline const vector& U() const;

        //- Return const access to the diameter
        inline scalar d() const;

        //- Return const access to the density
        inline scalar rho() const;

        //- Return const access to the mass flow rate
        inline scalar mDot() const;


    // Edit

        //- Return access to the position
        inline point& x();

        //- Return access to the velocity
        inline vector& U();

        //- Return access to the diameter
        inline scalar& d();

        //- Return access to the density
        inline scalar& rho();

        //- Return access to the mass flow rate
        inline scalar& mDot();


    // I-O

        // Friend Operators

            friend bool operator==
            (
                const streamParcelInjectionData& a,
                const streamParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const streamParcelInjectionData&, "
                        "const streamParcelInjectionData&"
                    ")"
                );

                return false;
            }

            friend bool operator!=
            (
                const streamParcelInjectionData& a,
                const streamParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const streamParcelInjectionData&, "
                        "const streamParcelInjectionData&"
                    ")"
                );

                return false;
            }

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const streamParcelInjectionData& data
            );

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                streamParcelInjectionData& data
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "streamParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
