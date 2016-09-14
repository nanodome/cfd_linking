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

\*---------------------------------------------------------------------------*/

#include "monteCarloParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::point& Foam::monteCarloParcelInjectionData::x() const
{
    return x_;
}


inline const Foam::vector& Foam::monteCarloParcelInjectionData::U() const
{
    return U_;
}


inline Foam::scalar Foam::monteCarloParcelInjectionData::d() const
{
    return d_;
}


inline Foam::scalar Foam::monteCarloParcelInjectionData::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::monteCarloParcelInjectionData::mDot() const
{
    return mDot_;
}


inline Foam::point& Foam::monteCarloParcelInjectionData::x()
{
    return x_;
}


inline Foam::vector& Foam::monteCarloParcelInjectionData::U()
{
    return U_;
}


inline Foam::scalar& Foam::monteCarloParcelInjectionData::d()
{
    return d_;
}


inline Foam::scalar& Foam::monteCarloParcelInjectionData::rho()
{
    return rho_;
}


inline Foam::scalar& Foam::monteCarloParcelInjectionData::mDot()
{
    return mDot_;
}


// ************************************************************************* //
