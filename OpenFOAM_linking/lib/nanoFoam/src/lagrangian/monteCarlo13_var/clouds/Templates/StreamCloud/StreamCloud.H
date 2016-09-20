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

Class
    Foam::StreamCloud

Description
    Templated base class for stream cloud

    - cloud function objects

    - particle forces, e.g.
      - buoyancy
      - drag
      - pressure gradient
      - ...

    - sub-models:
      - injection model
      - patch interaction model

SourceFiles
    StreamCloudI.H
    StreamCloud.C

\*---------------------------------------------------------------------------*/

#ifndef StreamCloud_H
#define StreamCloud_H

#include "streamCloud.H"

#include "XMLSampling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class StreamCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StreamCloud
:
    public CloudType,
    public streamCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef StreamCloud<CloudType> streamCloudType;


        typedef XMLSampling<StreamCloud<CloudType> > sampleType;

private:

    // Private data

        //- Cloud copy pointer
        autoPtr<StreamCloud<CloudType> > cloudCopyPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        StreamCloud(const StreamCloud&);

        //- Disallow default bitwise assignment
        void operator=(const StreamCloud&);


protected:

    // Protected data
	//- Hard coded sampling object for StreamCloud object
        sampleType xmlSampler_;

        // Initialisation

            //- Set cloud sub-models
            void setModels();

public:

    // Constructors

        //- Construct given carrier gas fields
        StreamCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
	    const volScalarField& T,
   	    const volScalarField& concSource,
            const volVectorField& gradT,	    
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        StreamCloud
        (
            StreamCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        StreamCloud
        (
            const fvMesh& mesh,
            const word& name,
            const StreamCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new StreamCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new StreamCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~StreamCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const StreamCloud& cloudCopy() const;

	    //- Access to xmlSampler methods
	    inline sampleType& xmlSampler();

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Evolve the cloud
            void evolve();

        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "StreamCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "StreamCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
