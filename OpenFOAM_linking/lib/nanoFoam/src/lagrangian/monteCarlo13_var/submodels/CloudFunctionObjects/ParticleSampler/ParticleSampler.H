/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::ParticleSampler

Description
    Function object to collect the parcel mass- and mass flow rate over a
    set of polygons.  The polygons are defined as lists of points.  If a
    parcel is 'collected', it is subsequently flagged to be removed from the
    domain.

SourceFiles
    ParticleSampler.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleSampler_H
#define ParticleSampler_H

#include "CloudFunctionObject.H"
#include "cylindricalCS.H"
#include "face.H"
#include "Switch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleSampler Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleSampler
:
    public CloudFunctionObject<CloudType>
{
public:

    enum modeType
    {
        mtPolygon,
        mtConcentricCircle,
        mtUnknown
    };


private:

    // Private Data

        // Typedefs

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Collector mode type
        modeType mode_;

        //- Index of parcel types to collect (-1 by default = all particles)
        const label parcelType_;

        //- Flag to remove collected particles
        Switch removeCollected_;

        //- List of points
        Field<point> points_;

        //- List of faces
        List<face> faces_;


        // Polygon collector

            //- Triangulation of faces
            List<List<face> > faceTris_;

        // Concentric circles collector

            //- Number of sectors per circle
            label nSector_;

            //- List of radii
            List<scalar> radius_;

            //- Cylindrical co-ordinate system
            cylindricalCS coordSys_;


        //- Face areas
        Field<scalar> area_;

        //- Polygon normal vector
        vector normal_;

        //- Remove mass of parcel travelling in opposite direction to normal_
        bool negateParcelsOppositeNormal_;

        //- Surface output format
        const word surfaceFormat_;

        //- Flag to indicate whether data should be reset/cleared on writing
        Switch resetOnWrite_;

        //- Total time
        scalar totalTime_;

        //- Mass storage
        List<scalar> mass_;

        //- Mass total storage
        List<scalar> massTotal_;

        //- Mass flow rate storage
        List<scalar> massFlowRate_;

        //- Flag to indicate whether data should be written to file
        Switch log_;

        //- Output file pointer
        autoPtr<OFstream> outputFilePtr_;

        //- Output directory
        fileName outputDir_;

        //- Last calculation time
        scalar timeOld_;


    // Private Member Functions

        //- Helper function to create log files
        void makeLogFile
        (
            const faceList& faces,
            const Field<point>& points,
            const Field<scalar>& area
        );

        //- Initialise polygon collectors
        void initPolygons();

        //- Initialise concentric circle collectors
        void initConcentricCircles();

        //- Collect parcels in polygon collectors
        label collectParcelPolygon
        (
            const point& position,
            const vector& U
        ) const;

        //- Collect parcels in concentric circle collectors
        label collectParcelConcentricCircles
        (
            const point& position,
            const vector& U
        ) const;

	// ----------- New Part ------------- //
	
	//- List of time for each data record
	List<DynamicList<scalar> > times_;

	//- List of output data of faces
	List<DynamicList<string> > faceData_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("particleSampler");


    // Constructors

        //- Construct from dictionary
        ParticleSampler(const dictionary& dict, CloudType& owner);

        //- Construct copy
        ParticleSampler(const ParticleSampler<CloudType>& pc);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new ParticleSampler<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleSampler();


    // Member Functions

        // Access

            //- Return const access to the reset on write flag
            inline const Switch& resetOnWrite() const;


        // Evaluation

            //- Post-move hook
            virtual void postMove
            (
                const parcelType& p,
                const label cellI,
                const scalar dt,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleSamplerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ParticleSampler.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
