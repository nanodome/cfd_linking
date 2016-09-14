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

#include "StreamParcel.H"
//#include "forceSuSp.H"
//#include "IntegrationScheme.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::StreamParcel<ParcelType>::maxTrackAttempts = 1;


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::StreamParcel<ParcelType>::setCellValues
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
                "void Foam::StreamParcel<ParcelType>::setCellValues"
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
}


template<class ParcelType>
template<class TrackData>
void Foam::StreamParcel<ParcelType>::cellValueSourceCorrection
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
void Foam::StreamParcel<ParcelType>::calc
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

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ = calcVelocity(td, dt, cellI, Re, muc_, mass0, Su, dUTrans, Spu);
   
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
const Foam::vector Foam::StreamParcel<ParcelType>::calcVelocity
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
    //Total velocity
    vector Unew = Uc_;    //Just take the velocity which is interpolated from the cell

    return Unew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::StreamParcel<ParcelType>::StreamParcel
(
    const StreamParcel<ParcelType>& p
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
    Tc_(p.Tc_),
    concSourcec_(p.concSourcec_),
    gradTc_(p.gradTc_)
{}


template<class ParcelType>
Foam::StreamParcel<ParcelType>::StreamParcel
(
    const StreamParcel<ParcelType>& p,
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
    Tc_(p.Tc_),
    concSourcec_(p.concSourcec_),
    gradTc_(p.gradTc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::StreamParcel<ParcelType>::move
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
	// Insert the old 2.2.x method call without start
        td.cloud().functions().postMove(p, cellI, dt, td.keepParticle);
    }

    const label cellI = p.cell();

    td.cloud().functions().nanoDomeCall(p, cellI, trackTime, td.keepParticle);
    //- Sampling call for XML sampling for this streamline solver
    td.cloud().xmlSampler().samplingCall(p, cellI, trackTime, td.keepParticle);

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
void Foam::StreamParcel<ParcelType>::hitFace(TrackData& td)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.cloud().functions().postFace(p, p.face(), td.keepParticle);
}


template<class ParcelType>
void Foam::StreamParcel<ParcelType>::hitFace(int& td)
{}


template<class ParcelType>
template<class TrackData>
bool Foam::StreamParcel<ParcelType>::hitPatch
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


template<class ParcelType>
template<class TrackData>
void Foam::StreamParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackData>
void Foam::StreamParcel<ParcelType>::hitWallPatch
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
void Foam::StreamParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParcelType>
void Foam::StreamParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    U_ = transform(T, U_);

    f_ = transform(T, f_);

    angularMomentum_ = transform(T, angularMomentum_);

    torque_ = transform(T, torque_);
}


template<class ParcelType>
void Foam::StreamParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


template<class ParcelType>
Foam::scalar Foam::StreamParcel<ParcelType>::wallImpactDistance
(
    const vector&
) const
{
    return 0.5*d_;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "StreamParcelIO.C"

// ************************************************************************* //
