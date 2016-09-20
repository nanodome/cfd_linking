/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "StreamCloud.H"
#include "IntegrationScheme.H"
#include "interpolation.H"
#include "subCycleTime.H"

//#include "InjectionModelList.H"
//#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::StreamCloud<CloudType>::setModels()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StreamCloud<CloudType>::StreamCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const volScalarField& T,
    const volScalarField& concSource,
    const volVectorField& gradT,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(cloudName, rho, U, mu, T, concSource, gradT, g, false),
    streamCloud(),
    xmlSampler_
    (
        *this,
        this->particleProperties_.subOrEmptyDict("sampling")
    )
{
}


template<class CloudType>
Foam::StreamCloud<CloudType>::StreamCloud
(
    StreamCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    streamCloud(),
    xmlSampler_(c.xmlSampler_)
{}


template<class CloudType>
Foam::StreamCloud<CloudType>::StreamCloud
(
    const fvMesh& mesh,
    const word& name,
    const StreamCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    streamCloud(),
    xmlSampler_(c.xmlSampler_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StreamCloud<CloudType>::~StreamCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CloudType>
void Foam::StreamCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<StreamCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::StreamCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::StreamCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::StreamCloud<CloudType>::evolve()
{
    if (this->solution_.canEvolve())
    {
        typename parcelType::template
            TrackingData<StreamCloud<CloudType> > td(*this);

        this->solve(td);
    }
}


template<class CloudType>
void Foam::StreamCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    typedef typename particle::TrackingData<StreamCloud<CloudType> > tdType;

    tdType td(*this);

    Cloud<parcelType>::template autoMap<tdType>(td, mapper);

    this->updateMesh();
}


// ************************************************************************* //
