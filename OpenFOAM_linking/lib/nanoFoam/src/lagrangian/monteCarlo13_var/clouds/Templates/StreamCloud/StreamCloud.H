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
// 
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

//#include "MonteCarloCloud.H"
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

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Dictionary of output properties
        IOdictionary outputProperties_;

        //- Solution properties
        cloudSolution solution_;

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;

        //- Sub-models dictionary
        const dictionary subModelProperties_;

        //- Random number generator - used by some injection routines
        cachedRandom rndGen_;

        //- Cell length scale
        scalarField cellLengthScale_;


        // References to the carrier gas fields

            //- Density [kg/m3]
            const volScalarField& rho_;

            //- Velocity [m/s]
            const volVectorField& U_;

            //- Dynamic viscosity [Pa.s]
            const volScalarField& mu_;

            //const volScalarField& Mf_;

	    //- added by Patrick
	    const volScalarField& T_;
	    const volScalarField& concSource_;
	    const volVectorField& gradT_;

        // Environmental properties
	    
            //- Gravity
            const dimensionedVector& g_;

            //- Averaged ambient domain pressure
            scalar pAmbient_;


	//- Hard coded sampling object for StreamCloud object
        sampleType xmlSampler_;


        // References to the cloud sub-models


        // Sources

            //- Momentum
            autoPtr<DimensionedField<vector, volMesh> > UTrans_;
           //- Temperatur  
            autoPtr<DimensionedField<scalar, volMesh> > TTrans_;
           //- Partikelanzahl	    
            autoPtr<DimensionedField<scalar, volMesh> > nParticleTrans_;	    
            autoPtr<DimensionedField<scalar, volMesh> > rhopTrans_;

  


            //- Coefficient for carrier phase U equation
            autoPtr<DimensionedField<scalar, volMesh> > UCoeff_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions


            //- Evolve the cloud
            template<class TrackData>
            void evolveCloud(TrackData& td);

            //- Reset state of cloud
            void cloudReset(StreamCloud<CloudType>& c);


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


            // References to the mesh and databases

                //- Return reference to the mesh
                inline const fvMesh& mesh() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return output properties dictionary
                inline const IOdictionary& outputProperties() const;

                //- Return non-const access to the output properties dictionary
                inline IOdictionary& outputProperties();

                //- Return const access to the solution properties
                inline const cloudSolution& solution() const;

                //- Return access to the solution properties
                inline cloudSolution& solution();

                //- Return the constant properties
                inline const typename parcelType::constantProperties&
                    constProps() const;

                    //- Return reference to the sub-models dictionary
                inline const dictionary& subModelProperties() const;


            // Cloud data

                //- Return reference to the random object
                inline cachedRandom& rndGen();

                //- Return the cell length scale
                inline const scalarField& cellLengthScale() const;


            // References to the carrier gas fields

                //- Return carrier gas velocity
                inline const volVectorField& U() const;

                //- Return carrier gas density
                inline const volScalarField& rho() const;

                //- Return carrier gas dynamic viscosity
                inline const volScalarField& mu() const;
	
		//- added by Patrick
		inline const volScalarField& T() const;
		inline const volScalarField& concSource() const;
		inline const volVectorField& gradT() const;
            // Environmental properties

                //- Gravity
                inline const dimensionedVector& g() const;

                //- Return const-access to the ambient pressure
                inline scalar pAmbient() const;

                //- Return reference to the ambient pressure
                inline scalar& pAmbient();


	    //- Access to xmlSampler methods
	    inline sampleType& xmlSampler();


            // Sources

                // Momentum

                    //- Return reference to momentum source
                    inline DimensionedField<vector, volMesh>& UTrans();

                    //- Return const reference to momentum source
                    inline const DimensionedField<vector, volMesh>&
                        UTrans() const;
			
                    inline DimensionedField<scalar, volMesh>& TTrans();

                
                    inline const DimensionedField<scalar, volMesh>&
                        TTrans() const;
                    
                    inline DimensionedField<scalar, volMesh>& nParticleTrans();

                
                    inline const DimensionedField<scalar, volMesh>&
                        nParticleTrans() const;	
                    
                    inline DimensionedField<scalar, volMesh>& rhopTrans();

                
                    inline const DimensionedField<scalar, volMesh>&
                        rhopTrans() const;	
			
                    
                    	
			
			
                   
			
		    
		   
         
                    			
                     //- Return coefficient for carrier phase U equation
                    inline DimensionedField<scalar, volMesh>& UCoeff();

                    //- Return const coefficient for carrier phase U equation
                    inline const DimensionedField<scalar, volMesh>&
                        UCoeff() const;

                    //- Return tmp momentum source term
                    inline tmp<fvVectorMatrix> SU(volVectorField& U) const;

                    // Temperatur von Partikeln auf Eulerfeld
                    inline const tmp<volScalarField> ST() const;
                    // Anzalh von Partikeln auf Eulerfeld
                    inline const tmp<volScalarField> SN() const;		    
                    // Dichte von Partikeln auf Eulerfeld		    
		    inline const tmp<volScalarField> Srhop() const; 
		    


                    
		    
        // Check

            //- Total number of parcels
            inline label nParcels() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total rotational kinetic energy in the system
            inline scalar rotationalKineticEnergyOfSystem() const;

            //- Penetration for fraction [0-1] of the current total mass
            inline scalar penetration(const scalar fraction) const;

            //- Mean diameter Dij
            inline scalar Dij(const label i, const label j) const;

            //- Max diameter
            inline scalar Dmax() const;


            // Fields

                //- Return the particle volume fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> theta() const;

                //- Return the particle mass fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> alpha() const;

                //- Return the particle effective density field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> rhoEff() const;


        // Cloud evolution functions


            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

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
