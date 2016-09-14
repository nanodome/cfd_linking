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

\*---------------------------------------------------------------------------*/

#include "MonteCarloParcel.H"
#include "forceSuSp.H"
#include "IntegrationScheme.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::MonteCarloParcel<ParcelType>::maxTrackAttempts = 1;


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    tetIndices tetIs = this->currentTetIndices();

    rhoc_ = td.rhoInterp().interpolate(this->position(), tetIs);

    if (rhoc_ < td.cloud().constProps().rhoMin())
    {
        if (debug)
        {
            WarningIn
            (
                "void Foam::MonteCarloParcel<ParcelType>::setCellValues"
                "("
                    "TrackData&, "
                    "const scalar, "
                    "const label"
                ")"
            )   << "Limiting observed density in cell " << cellI << " to "
                << td.cloud().constProps().rhoMin() <<  nl << endl;
        }

        rhoc_ = td.cloud().constProps().rhoMin();
    }

    Uc_ = td.UInterp().interpolate(this->position(), tetIs);

    muc_ = td.muInterp().interpolate(this->position(), tetIs);
    
    //Mfc_ = td.MfInterp().interpolate(this->position(), tetIs);

    //epsilonc_ = td.epsilonInterp().interpolate(this->position(), tetIs);

    //kc_ = td.kInterp().interpolate(this->position(), tetIs);

    //YTTIPc_ = td.YTTIPInterp().interpolate(this->position(), tetIs);  
      
    //- added by Patrick

    Tc_ = td.TInterp().interpolate(this->position(), tetIs); 

    concSourcec_ = td.concSourceInterp().interpolate(this->position(), tetIs);

    gradTc_ = td.gradTInterp().interpolate(this->position(), tetIs); 

        // Apply dispersion components to carrier phase velocity
    Uc_ = td.cloud().dispersion().update
    (
        dt,
        cellI,
        U_,
        Uc_,
        UTurb_,
        tTurb_
    );
}



template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::calcPBE
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
   scalar kt = 1.0;	//sectional model value would be 1.0e-03
   label nsteps = 0;
   scalar parDt; //= dt/nsteps;
   scalar dtTot = 0.0;
   //scalar Cun_ = 0.0;

   //for (dti=1; dti <= nsteps; dti++)

   do {
	   parDt =  max( kt / max(beta_, 1.0e-40) / max(N_, 1.0e-40) , 1e-8 );

	   parDt = min ( parDt, dt - dtTot);	//shortenig the last timestep to match dt
	   dtTot += parDt;
	   nsteps += 1;

   //- modified by Patrick
           
	   //Info<<"+++++++++++++++++++++++++++++++++++++++++++++++ "<<endl;       
          
	   N_ += concSourcec_*6.022e23*parDt -parDt*0.5*beta_*sqr(N_);     
           N_ = max(N_, 0.0); 
    
           V_ += concSourcec_*6.022e23*parDt*6.544985e-29;		//3.351e-29;
           V_ = max(V_, 0.0); 
	    
	   //area density of fused, spherical particles
           if (N_ > 0.0)
           {
               a_ = Foam::pow(V_/(6.544985e-29*N_),2.0/3.0) * N_ * 7.8539816e-19;	//5.026548e-19;
	   }
           a_ = max(a_, 0.0);  

           //- modified by Patrick
	   //- in this try there is no sintering  
           
	   //if(tau_ > 0.0 && tau_< 10000 ) //>= 0.001*A_)
           //{  
           //     A_ +=RRp_*dt*5.026548e-19-(A_-N_*a_)*dt/tau_;  
	   //    
           //}
           //else
           //{  
                //A_ +=concSourcec_*6.022e23*parDt*7.8539816e-19;		//5.026548e-19;
		A_ +=concSourcec_*6.022e23*parDt*7.8539816e-19 -(A_ - a_)/tau_;
	   //    
           //}	    
	      
           A_ = max(A_, a_);   
     
            if (N_ > 0.0 ) 
            {
                v_ = V_ / N_;
            }
           //v_ = max(v_, 0.0); 
	
            if (N_ > 0.0) 
            {	    
	          //dp_ = Foam::pow(6.0/3.1415926*v_, 1.0/3.0);   //diameter of an aggregate
		  dp_ = 6.0 * V_ / A_;
	    }
           //dp_ = max(dp_, 0.0);

	   //collision radius of an aggregate
           if (A_ > 0.0 && V_ > 0.0 ) 
           {
               //rc_ = 0.5 * (6.0 * V_ / A_) * Foam::pow(V_ / N_ * 6.0 / (3.1415926 * Foam::pow(6.0 * V_ / A_, 3.0)), 1.0/1.8);
	
               rc_ = 3.0 * V_ / A_ * Foam::pow(Foam::pow(A_, 3.0) / 36.0 / 3.1415926 /(Foam::sqr(V_ ) * N_), 1.0/1.8);
 
	   // 0.5 - radius    primary particle size      number of primary particles per aggregate
	   }
           //rc_ = max(rc_, 0.0); 

           //particle diffusivity
           if (rc_ > 0.0) 
           {
              //Kn_ = muc_ / 1.0e05 * Foam::sqrt(3.1415926 * Tc_ * 8.3144 / (2.0 * 0.032)) / rc_;
              Kn_ = muc_ / 1.0e05 * Foam::sqrt(3.1415926 * Tc_ * 1.3806485e-23 * 6.022e23 / (2.0 * 0.028)) / rc_;
	      Cun_ = (5.0 + 4.0 * Kn_ + 6.0 * Foam::sqr(Kn_) + 18.0 * Foam::pow(Kn_,3.0)) / (5.0 - Kn_ + (8.0 + 3.1415926) * Foam::sqr(Kn_));
	      Dp_= 1.3806504e-23*Tc_ / 6.0 / 3.1415926 / muc_ /rc_; //* Cun_;
           }
            //Dp_ = max(Dp_, 0.0); 

           if (v_ > 0.0)
           {
               cn_=Foam::pow(1.3806504e-23*8.0*Tc_/(3.1415926*rho_*v_), 1.0/2.0);   
	    
           }
           //cn_ = max(cn_, 0.0);

           if (cn_ > 0.0)
           {
               l_=8.0*Dp_/3.1415926/cn_;
	    
           }
           //l_ = max(l_, 0.0);

           if (l_ > 0.0 && rc_ > 0.0) 
          {
      	       gp_ = 1.0/(6.0*rc_*l_)*(Foam::pow((2.0*rc_+l_), 3.0)-Foam::pow(4.0*rc_*rc_+l_*l_, 3.0/2.0))-2.0*rc_;
	
          }
          //gp_ = max(gp_, 0.0); 
   
          if (cn_ > 0.0 && rc_ > 0.0) 
          {
	     beta_ = 8.0*3.1415926*Dp_*rc_*Foam::pow((rc_/(2.0*rc_+Foam::pow(2.0, 1.0/2.0)*gp_)+Foam::pow(2.0, 1.0/2)*Dp_/(cn_*rc_)), -1.0);

          //beta_ = max(beta_, 0.0);
	  }

          if (A_ > 0.0 ) 
          {
	      //tau_ = 6.5e-17 * 6.0 * V_ / A_ * Foam::exp(83000/(8.3144 * max(650.0, Tc_)));	//Titanium oxid 
	      tau_ = 8.0e16 * Tc_ * Foam::pow(6.0 * V_ / A_, 4.0) * Foam::exp(30000.0 / Tc_);	//Iron oxid 
          }
          //tau_ = max(tau_, 0.0); 
   
 
  } while(dtTot < dt || nsteps >= 1000);

  d_ = Foam::pow(6.0/3.1415926*v_, 1.0/3.0);	//For drag and wight forces 

  if (nsteps == 1000)
  {
      Info << "Warning: maximum step count for PBE reached. Particle " << this->origId() << endl;
  }
}


template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    Uc_ += td.cloud().UTrans()[cellI]/massCell(cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();

    // Reynolds number
    const scalar Re = this->Re(U_, d_, rhoc_, muc_);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    // - modified by Patrick
    vector Su = vector::zero;
    // thermophoretic force on particle 
    //scalar FreePath = muc_ / 1.0e05 * Foam::sqrt(3.1415926 * Tc_ * 1.3806485e-23 * 6.022e23 / (2.0 * 0.028));
    //Su = -1.0e5 * FreePath * Foam::sqr(d_) * gradTc_ / Tc_;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ = calcVelocity(td, dt, cellI, Re, muc_, mass0, Su, dUTrans, Spu);	//U from forces
    //U_ = Uc_;		// U from velocity field
   
    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().solution().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Update momentum transfer coefficient
        td.cloud().UCoeff()[cellI] += np0*Spu;
    }
}


template<class ParcelType>
template<class TrackData>
const Foam::vector Foam::MonteCarloParcel<ParcelType>::calcVelocity
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar mu,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu
) const
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::forceType forceType;

    //- modified by Patrick

    // Momentum source due to particle forces
    const parcelType& p = static_cast<const parcelType&>(*this);

    //- modified by Patrick
    // Velocity calculation due to gas velovity and thermophoretic forces (in newest version, not here)
     
    //Random walk caused by brownian motion
    cachedRandom& rnd = td.cloud().rndGen();

    scalar sigma = 0.0;
    if (dt > 1.0e-10)
    {
        sigma = Foam::sqrt(2.0 * Dp_ / dt);
    }


    vector dir = 2.0*rnd.sample01<vector>() - vector::one;
    dir /= mag(dir) + SMALL;

    // Numerical Recipes... Ch. 7. Random Numbers...
    scalar x1 = 0.0;
    scalar x2 = 0.0;
    scalar rsq = 10.0;
    while ((rsq > 1.0) || (rsq == 0.0))
    {
        x1 = 2.0*rnd.sample01<scalar>() - 1.0;
        x2 = 2.0*rnd.sample01<scalar>() - 1.0;
        rsq = x1*x1 + x2*x2;
    }

    scalar fac = sqrt(-2.0*log(rsq)/rsq);

    fac *= mag(x1);

    vector UB = sigma * fac * dir;
    
    //Total velocity
    vector Unew = Uc_ +UB;    //Summ of resulting velocities of particle

    return Unew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::MonteCarloParcel<ParcelType>::MonteCarloParcel
(
    const MonteCarloParcel<ParcelType>& p
)
:
    ParcelType(p),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    f_(p.f_),
    angularMomentum_(p.angularMomentum_),
    torque_(p.torque_),
    rho_(p.rho_),
    age_(p.age_),
    cV_(p.cV_),
    //- modified by Patrick
    //Mfp_(p.Mfp_), 
    //Tp_(p.Tp_),
    //YTTIPp_(p.YTTIPp_), 
    //RRp_(p.RRp_),
    N_(p.N_),
    A_(p.A_),
    V_(p.V_),
    a_(p.a_),
    v_(p.v_),    
    dp_(p.dp_),
    rc_(p.rc_),
    cn_(p.cn_),
    Dp_(p.Dp_),
    Kn_(p.Kn_),
    Cun_(p.Cun_),
    l_(p.l_),
    gp_(p.gp_),
    beta_(p.beta_),
    tau_(p.tau_), 
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_),
    //Mfc_(p.Mfc_),
    //epsilonc_(p.epsilonc_),
    //kc_(p.kc_),
    //YTTIPc_(p.YTTIPc_)
    //- added by Patrick
    Tc_(p.Tc_),
    concSourcec_(p.concSourcec_),
    gradTc_(p.gradTc_)
{}


template<class ParcelType>
Foam::MonteCarloParcel<ParcelType>::MonteCarloParcel
(
    const MonteCarloParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    f_(p.f_),
    angularMomentum_(p.angularMomentum_),
    torque_(p.torque_),
    rho_(p.rho_),
    age_(p.age_),
    cV_(p.cV_),
    //- modiefied by Patrick
    //Mfp_(p.Mfp_),
    //Tp_(p.Tp_), 
    //YTTIPp_(p.YTTIPp_), 
    //RRp_(p.RRp_), 
    N_(p.N_),
    A_(p.A_),
    V_(p.V_),
    a_(p.a_),
    v_(p.v_),    
    dp_(p.dp_),
    rc_(p.rc_),
    cn_(p.cn_),
    Dp_(p.Dp_),
    Kn_(p.Kn_),
    Cun_(p.Cun_),
    l_(p.l_),
    gp_(p.gp_),
    beta_(p.beta_),
    tau_(p.tau_),     
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    rhoc_(p.rhoc_),
    Uc_(p.Uc_),
    muc_(p.muc_),
    //Mfc_(p.Mfc_),
    //epsilonc_(p.epsilonc_),
    //kc_(p.kc_),
    //YTTIPc_(p.YTTIPc_)
    //- added by Patrick
    Tc_(p.Tc_),
    concSourcec_(p.concSourcec_),
    gradTc_(p.gradTc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::MonteCarloParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
    const scalarField& cellLengthScale = td.cloud().cellLengthScale();
    const scalar maxCo = td.cloud().solution().maxCo();

    scalar tEnd = (1.0 - p.stepFraction())*trackTime;
    scalar dtMax = trackTime;
    if (td.cloud().solution().transient())
    {
        dtMax *= maxCo;
    }

    bool tracking = true;
    label nTrackingStalled = 0;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        // Apply correction to position for reduced-D cases
        meshTools::constrainToMeshCentre(mesh, p.position());

        const point start(p.position());

        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Cache the parcel current cell as this will change if a face is hit
        const label cellI = p.cell();

        const scalar magU = mag(U_);
        if (p.active() && tracking && (magU > ROOTVSMALL))
        {
            const scalar d = dt*magU;
            const scalar dCorr = min(d, maxCo*cellLengthScale[cellI]);
            dt *=
                dCorr/d
               *p.trackToFace(p.position() + dCorr*U_/magU, td);
        }

        tEnd -= dt;

        scalar newStepFraction = 1.0 - tEnd/trackTime;

        if (tracking)
        {
            if
            (
                mag(p.stepFraction() - newStepFraction)
              < particle::minStepFractionTol
            )
            {
                nTrackingStalled++;

                if (nTrackingStalled > maxTrackAttempts)
                {
                    tracking = false;
                }
            }
            else
            {
                nTrackingStalled = 0;
            }
        }

        p.stepFraction() = newStepFraction;

        bool calcParcel = true;
        if (!tracking && td.cloud().solution().steadyState())
        {
            calcParcel = false;
        }

        // Avoid problems with extremely small timesteps
        if ((dt > ROOTVSMALL) && calcParcel)
        {
            // Update cell based properties
            p.setCellValues(td, dt, cellI);

	    //- added by Patrick
	    p.calcPBE(td, dt, cellI);

            if (td.cloud().solution().cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(td, dt, cellI);
            }

            p.calc(td, dt, cellI);
        }

        if (p.onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
            {
                td.switchProcessor = true;
            }
        }

        p.age() += dt;

	// Out comment new method call with dt and start
        //td.cloud().functions().postMove(p, cellI, dt, start, td.keepParticle);
	// Insert the old method call without start
        td.cloud().functions().postMove(p, cellI, dt, td.keepParticle);
    }

    const label cellI = p.cell();

    td.cloud().functions().nanoDomeCall(p, cellI, trackTime, td.keepParticle);

    return td.keepParticle;
}


// old 2.2.x move method

//template<class ParcelType>
//template<class TrackData>
//bool Foam::MonteCarloParcel<ParcelType>::move
//(
//    TrackData& td,
//    const scalar trackTime
//)
//{
//    typename TrackData::cloudType::parcelType& p =
//        static_cast<typename TrackData::cloudType::parcelType&>(*this);
//
//    td.switchProcessor = false;
//    td.keepParticle = true;
//
//    const polyMesh& mesh = td.cloud().pMesh();
//    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
//    const scalarField& V = mesh.cellVolumes();
//    const scalar maxCo = td.cloud().solution().maxCo();
//
//    scalar tEnd = (1.0 - p.stepFraction())*trackTime;
//    const scalar dtMax = tEnd;
//
//    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
//    {
//	
//        // Apply correction to position for reduced-D cases
//        meshTools::constrainToMeshCentre(mesh, p.position());
//
//        // Set the Lagrangian time-step
//        scalar dt = min(dtMax, tEnd);
//
//	//Info << "Actual time step size:" << dt << endl;
//
//        // Remember which cell the parcel is in since this will change if
//        // a face is hit
//        const label cellI = p.cell();
//
//        const scalar magU = mag(U_);
//        if (p.active() && magU > ROOTVSMALL)
//        {
//            const scalar d = dt*magU;
//            const scalar dCorr = min(d, maxCo*cbrt(V[cellI]));
//            dt *=
//                dCorr/d
//               *p.trackToFace(p.position() + dCorr*U_/magU, td);
//
//	    //Info << "d: " << d << "  dCorr: " << dCorr << "  dt: " << dt << endl;
//        }
//	
//        tEnd -= dt;
//        p.stepFraction() = 1.0 - tEnd/trackTime;
//
//	//Info << "New tEnd: " << tEnd << "  stepFraction: " << p.stepFraction() << endl;
//
//        // Avoid problems with extremely small timesteps
//        if (dt > ROOTVSMALL)
//        {
//	    //Info << "Update values now. U before: " << U_ << "  pos: " << p.position() << endl;
//
//            // Update cell based properties
//            p.setCellValues(td, dt, cellI);
//
//	    //- added by Patrick
//	    p.calcPBE(td, dt, cellI);
//
//            if (td.cloud().solution().cellValueSourceCorrection())
//            {
//                p.cellValueSourceCorrection(td, dt, cellI);
//            }
//
//            p.calc(td, dt, cellI);
//
//	    //Info << "Update values complete U new: " << U_ << "  pos: " << p.position() << endl;
//	     
//        }
//
//        if (p.onBoundary() && td.keepParticle)
//        {
//            if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
//            {
//                td.switchProcessor = true;
//            }
//        }
//
//	//Info << "Switch prozessor? " << td.switchProcessor << endl;
//
//        p.age() += dt;
//
//        td.cloud().functions().postMove(p, cellI, dt, td.keepParticle);
//
//	//Info << "Post move done!" << endl;
//    }
//
//    //Info << "Move cycle finished!" << endl;
//
//    const label cellI = p.cell();
//    
//    td.cloud().functions().nanoDomeCall(p, cellI, trackTime, td.keepParticle);
//
//    //Info << "NanoDome call done!" << endl;
//
//    return td.keepParticle;
//}


template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::hitFace(TrackData& td)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.cloud().functions().postFace(p, p.face(), td.keepParticle);
}


template<class ParcelType>
void Foam::MonteCarloParcel<ParcelType>::hitFace(int& td)
{}


template<class ParcelType>
template<class TrackData>
bool Foam::MonteCarloParcel<ParcelType>::hitPatch
(
    const polyPatch& pp,
    TrackData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    // Invoke post-processing model
    td.cloud().functions().postPatch
    (
        p,
        pp,
        trackFraction,
        tetIs,
        td.keepParticle
    );

    // Invoke surface film model
    if (td.cloud().surfaceFilm().transferParcel(p, pp, td.keepParticle))
    {
        // All interactions done
        return true;
    }
    else
    {
        // Invoke patch interaction model
        return td.cloud().patchInteraction().correct
        (
            p,
            pp,
            td.keepParticle,
            trackFraction,
            tetIs
        );
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td,
    const tetIndices&
)
{
    // Wall interactions handled by generic hitPatch function
}


template<class ParcelType>
template<class TrackData>
void Foam::MonteCarloParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParcelType>
void Foam::MonteCarloParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    U_ = transform(T, U_);

    f_ = transform(T, f_);

    angularMomentum_ = transform(T, angularMomentum_);

    torque_ = transform(T, torque_);
}


template<class ParcelType>
void Foam::MonteCarloParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


template<class ParcelType>
Foam::scalar Foam::MonteCarloParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    return 0.5*d_;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "MonteCarloParcelIO.C"

// ************************************************************************* //
