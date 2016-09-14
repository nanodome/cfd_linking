/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::MonteCarloParcel

Description
    MonteCarlo parcel class with rotational motion (as spherical
    particles only) and one/two-way coupling with the continuous
    phase.

    Sub-models include:
    - drag
    - turbulent dispersion
    - wall interactions

SourceFiles
    MonteCarloParcelI.H
    MonteCarloParcel.C
    MonteCarloParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef MonteCarloParcel_H
#define MonteCarloParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolation.H"

// #include "ParticleForceList.H" // TODO

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class MonteCarloParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const MonteCarloParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class MonteCarloParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class MonteCarloParcel
:
    public ParcelType
{
public:

    //- Class to hold monteCarlo particle constant properties
    class constantProperties
    {
        // Private data

            //- Constant properties dictionary
            const dictionary dict_;

            //- Parcel type id - used for post-processing to flag the type
            //  of parcels issued by this cloud
            label parcelTypeId_;

            //- Minimum density [kg/m3]
            scalar rhoMin_;

            //- Particle density [kg/m3] (constant)
            scalar rho0_;

            //- Minimum particle mass [kg]
            scalar minParticleMass_;

            //- Young's modulus [N/m2]
            scalar youngsModulus_;

            //- Poisson's ratio
            scalar poissonsRatio_;

            scalar Cy_;	    


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Constructor from dictionary
            constantProperties
            (
                const dictionary& parentDict,
                const bool readFields = true
            );

            //- Construct from components
            constantProperties
            (
                const label parcelTypeId,
                const scalar rhoMin,
                const scalar rho0,
                const scalar minParticleMass,
                const scalar youngsModulus,
                const scalar poissonsRatio,
                const scalar Cy  
            );


        // Member functions

            //- Return const access to the constant properties dictionary
            inline const dictionary& dict() const;

            //- Return const access to the parcel type id
            inline label parcelTypeId() const;

            //- Return const access to the minimum density
            inline scalar rhoMin() const;

            //- Return const access to the particle density
            inline scalar rho0() const;

            //- Return const access to the minimum particle mass
            inline scalar minParticleMass() const;

            //- Return const access to Young's Modulus
            inline scalar youngsModulus() const;

            //- Return const access to Poisson's ratio
            inline scalar poissonsRatio() const;
	    
            inline scalar Cy() const;	  	    
    };


    template<class CloudType>
    class TrackingData
    :
        public ParcelType::template TrackingData<CloudType>
    {
    public:

        enum trackPart
        {
            tpVelocityHalfStep,
            tpLinearTrack,
            tpRotationalTrack
        };


    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Density interpolator
                autoPtr<interpolation<scalar> > rhoInterp_;

                //- Velocity interpolator
                autoPtr<interpolation<vector> > UInterp_;
		

                //- Dynamic viscosity interpolator
                autoPtr<interpolation<scalar> > muInterp_;
		//- modified by Patrick
                //autoPtr<interpolation<scalar> > MfInterp_;		
                //autoPtr<interpolation<scalar> > epsilonInterp_;		
                //autoPtr<interpolation<scalar> > kInterp_;
                //autoPtr<interpolation<scalar> > YTTIPInterp_;
		//- added by Patrick
		autoPtr<interpolation<scalar> > TInterp_;		
		autoPtr<interpolation<scalar> > concSourceInterp_;
		autoPtr<interpolation<vector> > gradTInterp_;

            //- Local gravitational or other body-force acceleration
            const vector& g_;

            // label specifying which part of the integration
            // algorithm is taking place
            trackPart part_;


    public:

        // Constructors

            //- Construct from components
            inline TrackingData
            (
                CloudType& cloud,
                trackPart part = tpLinearTrack
            );


        // Member functions

            //- Return conat access to the interpolator for continuous
            //  phase density field
            inline const interpolation<scalar>& rhoInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase velocity field
            inline const interpolation<vector>& UInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase dynamic viscosity field
            inline const interpolation<scalar>& muInterp() const;
            
	    //- modified by Patrick
            //inline const interpolation<scalar>& MfInterp() const;
            //inline const interpolation<scalar>& epsilonInterp() const;
	    //inline const interpolation<scalar>& kInterp() const;
	    //inline const interpolation<scalar>& YTTIPInterp() const;
	    //-added by Patrick 
	    inline const interpolation<scalar>& TInterp() const;	    
    	    inline const interpolation<scalar>& concSourceInterp() const;
	    inline const interpolation<vector>& gradTInterp() const; 
	    
            // Return const access to the gravitational acceleration vector
            inline const vector& g() const;

            //- Return the part of the tracking operation taking place
            inline trackPart part() const;

            //- Return access to the part of the tracking operation taking place
            inline trackPart& part();
    };


    //- Number of particle tracking attempts before we assume that it stalls
    static label maxTrackAttempts;


protected:
    
        // Parcel properties

            //- Active flag - tracking inactive when active = false
            bool active_;

            //- Parcel type id
            label typeId_;

            //- Number of particles in Parcel
            scalar nParticle_;

            //- Diameter [m]
            scalar d_;

            //- Target diameter [m]
            scalar dTarget_;

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Force on particle due to collisions [N]
            vector f_;

            //- Angular momentum of Parcel in global reference frame
            // [kg m2/s]
            vector angularMomentum_;

            //- Torque on particle due to collisions in global
            //  reference frame [Nm]
            vector torque_;

            //- Density [kg/m3]
            scalar rho_;

            //- Age [s]
            scalar age_;
	    scalar cV_;
	    //scalar Mfp_;
	    
            //scalar Tp_;	
            //scalar YTTIPp_;	    
            //scalar RRp_;	    		    
            scalar N_;
            scalar A_;	    
            scalar V_;	    
            scalar a_;	    
            scalar v_;
            scalar dp_;	    
            scalar rc_;	    
            scalar cn_;
	    scalar Dp_;
	    scalar Kn_;
            scalar Cun_;	    
            scalar l_;	    
            scalar gp_;	    
            scalar beta_;	    
            scalar tau_;	    
	    
	    
            //- Time spent in turbulent eddy [s]
            scalar tTurb_;

            //- Turbulent velocity fluctuation [m/s]
            vector UTurb_;


        // Cell-based quantities

            //- Density [kg/m3]
            scalar rhoc_;

            //- Velocity [m/s]
            vector Uc_;

            //- Viscosity [Pa.s]
            scalar muc_;
	
	    //- modified by Patrick
            //scalar Mfc_;
            //scalar epsilonc_;    

            //scalar kc_; 
            //scalar YTTIPc_;
	    scalar Tc_;	 
	    //- added by Patrick
	    scalar concSourcec_;
	    vector gradTc_;
	    
	    
    // Protected Member Functions

        //- Calculate new particle velocity
        template<class TrackData>
        const vector calcVelocity
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar mu,           // local carrier viscosity
            const scalar mass,         // mass
            const vector& Su,          // explicit particle momentum source
            vector& dUTrans,           // momentum transfer to carrier
            scalar& Spu                // linearised drag coefficient
        ) const;


public:

    // Static data members

        //- Runtime type information
        TypeName("MonteCarloParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " active"
          + " typeId"
          + " nParticle"
          + " d"
          + " dTarget "
          + " (Ux Uy Uz)"
          + " (fx fy fz)"
          + " (angularMomentumx angularMomentumy angularMomentumz)"
          + " (torquex torquey torquez)"
          + " rho"
          + " age"
	  + " cV"
          //+ " Mfp"
          //+ " Tp"
          //+ " YTTIPp"
          //+ " RRp"
          + " N"	  
          + " A"	  
          + " V"	  
          + " a"	  
          + " v"	  
          + " dp"	  
          + " rc"	
          + " cn"	
	  + " Dp"	
	  + " Kn"	
	  + " Cun"	
          + " l"	
          + " gp"	
          + " beta"	
          + " tau"	  
          + " tTurb"
          + " (UTurbx UTurby UTurbz)"
	  + " Tc"
	  + " concSourcec"
	  + " (gradTcx gradTcy gradTcz)"
        );


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline MonteCarloParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline MonteCarloParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            //const scalar Mfp0,	
            //const scalar YTTIPp0,		 
            const constantProperties& constProps
        );

        //- Construct from Istream
        MonteCarloParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        MonteCarloParcel(const MonteCarloParcel& p);

        //- Construct as a copy
        MonteCarloParcel(const MonteCarloParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new MonteCarloParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new MonteCarloParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<MonteCarloParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<MonteCarloParcel<ParcelType> >
                (
                    new MonteCarloParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to active flag
            inline bool active() const;

            //- Return const access to type id
            inline label typeId() const;

            //- Return const access to number of particles
            inline scalar nParticle() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to target diameter
            inline scalar dTarget() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to force
            inline const vector& f() const;

            //- Return const access to angular momentum
            inline const vector& angularMomentum() const;

            //- Return const access to torque
            inline const vector& torque() const;

            //- Return const access to density
            inline scalar rho() const;

            //- Return const access to the age
            inline scalar age() const;
            inline scalar cV() const;
	    //- mpdified by Patrick
            //inline scalar Mfp() const; 

            //inline scalar Tp() const; 

            //inline scalar YTTIPp() const; 	    
            //inline scalar RRp() const; 
            inline scalar N() const; 	    
            inline scalar A() const; 	    
            inline scalar V() const; 	    
            inline scalar a() const; 	    
            inline scalar v() const; 	    
            inline scalar dp() const; 	    
            inline scalar rc() const; 	    
            inline scalar cn() const;
            inline scalar Dp() const;
	    inline scalar Kn() const;
	    inline scalar Cun() const;
            inline scalar l() const; 	    
            inline scalar gp() const; 	    
            inline scalar beta() const; 	    
            inline scalar tau() const; 	    
            //- Return const access to time spent in turbulent eddy
            inline scalar tTurb() const;

            //- Return const access to turbulent velocity fluctuation
            inline const vector& UTurb() const;

            //- Return const access to carrier density [kg/m3]
            inline scalar rhoc() const;

            //- Return const access to carrier velocity [m/s]
            inline const vector& Uc() const;

            //- Return const access to carrier viscosity [Pa.s]
            inline scalar muc() const;
	   
	    //- modified by Patrick 
            //inline scalar Mfc() const;

            //inline scalar epsilonc() const;
	   
            //inline scalar kc() const;

            //inline scalar YTTIPc() const;
	    
	    //- added by Patrick

            inline scalar Tc() const;
            inline scalar concSourcec() const;
	    inline const vector& gradTc() const;	    
	    
	    
        // Edit

            //- Return const access to active flag
            inline bool& active();

            //- Return access to type id
            inline label& typeId();

            //- Return access to number of particles
            inline scalar& nParticle();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to target diameter
            inline scalar& dTarget();

            //- Return access to velocity
            inline vector& U();

            //- Return access to force
            inline vector& f();

            //- Return access to angular momentum
            inline vector& angularMomentum();

            //- Return access to torque
            inline vector& torque();

            //- Return access to density
            inline scalar& rho();

            //- Return access to the age
            inline scalar& age();
	    inline scalar& cV();
	    //- modified by Patrick

            //inline scalar& Mfp();
            //inline scalar& Tp();
            //inline scalar& YTTIPp();	    
            //inline scalar& RRp();
            inline scalar& N();	    
            inline scalar& A();	    
            inline scalar& V();	    
            inline scalar& a();	    
            inline scalar& v();	    
            inline scalar& dp();	    
            inline scalar& rc();	    
            inline scalar& cn();
            inline scalar& Dp();
	    inline scalar& Kn();
	    inline scalar& Cun();
            inline scalar& l();	    
            inline scalar& gp();	    
            inline scalar& beta();	    
            inline scalar& tau();	    
            //- Return access to time spent in turbulent eddy
            inline scalar& tTurb();

            //- Return access to turbulent velocity fluctuation
            inline vector& UTurb();


        // Helper functions

            //- Return the index of the face to be used in the interpolation
            //  routine
            inline label faceInterpolation() const;

            //- Cell owner mass
            inline scalar massCell(const label cellI) const;

            //- Particle mass
            inline scalar mass() const;

            //- Particle moment of inertia around diameter axis
            inline scalar momentOfInertia() const;

            //- Particle angular velocity
            inline vector omega() const;

            //- Particle volume
            inline scalar volume() const;

            //- Particle volume for a given diameter
            inline static scalar volume(const scalar d);

            //- Particle projected area
            inline scalar areaP() const;

            //- Projected area for given diameter
            inline static scalar areaP(const scalar d);

            //- Particle surface area
            inline scalar areaS() const;

            //- Surface area for given diameter
            inline static scalar areaS(const scalar d);

            //- Reynolds number
            inline scalar Re
            (
                const vector& U,        // particle velocity
                const scalar d,         // particle diameter
                const scalar rhoc,      // carrier density
                const scalar muc        // carrier dynamic viscosity
            ) const;

            //- Weber number
            inline scalar We
            (
                const vector& U,        // particle velocity
                const scalar d,         // particle diameter
                const scalar rhoc,      // carrier density
                const scalar sigma      // particle surface tension
            ) const;

            //- Eotvos number
            inline scalar Eo
            (
                const vector& a,        // acceleration
                const scalar d,         // particle diameter
                const scalar sigma      // particle surface tension
            ) const;


        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );
  
	    
	    
            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );
            template<class TrackData>
            void calcPBE
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );


        // Tracking

            //- Move the parcel
            template<class TrackData>
            bool move(TrackData& td, const scalar trackTime);


        // Patch interactions

            //- Overridable function to handle the particle hitting a face
            //  without trackData
            void hitFace(int& td);

            //- Overridable function to handle the particle hitting a face
            template<class TrackData>
            void hitFace(TrackData& td);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackData>
            bool hitPatch
            (
                const polyPatch& p,
                TrackData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackData>
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackData>
            void hitWallPatch
            (
                const wallPolyPatch&,
                TrackData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            template<class TrackData>
            void hitPatch
            (
                const polyPatch&,
                TrackData& td
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);

            //- The nearest distance to a wall that the particle can be
            //  in the n direction
            virtual scalar wallImpactDistance(const vector& n) const;


        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const MonteCarloParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MonteCarloParcelI.H"
#include "MonteCarloParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MonteCarloParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
