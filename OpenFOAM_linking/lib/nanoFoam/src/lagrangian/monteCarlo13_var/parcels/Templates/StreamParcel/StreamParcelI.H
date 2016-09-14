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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline
Foam::StreamParcel<ParcelType>::constantProperties::constantProperties()
:
    dict_(dictionary::null),
    parcelTypeId_(-1),
    rhoMin_(0.0),
    rho0_(0.0),
    minParticleMass_(0.0),
    youngsModulus_(0.0),
    poissonsRatio_(0.0),
    Cy_(0.0)  
{}


template<class ParcelType>
inline Foam::StreamParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    dict_(cp.dict_),
    parcelTypeId_(cp.parcelTypeId_),
    rhoMin_(cp.rhoMin_),
    rho0_(cp.rho0_),
    minParticleMass_(cp.minParticleMass_),
    youngsModulus_(cp.youngsModulus_),
    poissonsRatio_(cp.poissonsRatio_),
    Cy_(cp.Cy_)
{}


template<class ParcelType>
inline Foam::StreamParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict,
    const bool readFields
)
:
    dict_(parentDict.subOrEmptyDict("constantProperties")),
    parcelTypeId_(1),
    rhoMin_(1e-15),
    rho0_(0.0),
    minParticleMass_(1e-15),
    youngsModulus_(0.0),
    poissonsRatio_(0.0),
    Cy_(0.0) 
{
    if (readFields)
    {
        if (dict_.readIfPresent("parcelTypeId", parcelTypeId_))
        {
            Info<< "    employing parcel parcelTypeId of " << parcelTypeId_
                << endl;
        }
        if (dict_.readIfPresent("rhoMin", rhoMin_))
        {
            Info<< "    employing parcel rhoMin of " << rhoMin_ << endl;
        }
        if (dict_.readIfPresent("minParticleMass", minParticleMass_))
        {
            Info<< "    employing parcel minParticleMass of "
                << minParticleMass_ << endl;
        }

        dict_.lookup("rho0") >> rho0_;
        dict_.lookup("youngsModulus") >> youngsModulus_;
        dict_.lookup("poissonsRatio") >> poissonsRatio_;
        dict_.lookup("Cy") >> Cy_;	
    }
}

template<class ParcelType>
inline Foam::StreamParcel<ParcelType>::constantProperties::constantProperties
(
    const label parcelTypeId,
    const scalar rhoMin,
    const scalar rho0,
    const scalar minParticleMass,
    const scalar youngsModulus,
    const scalar poissonsRatio,
    const scalar Cy 
)
:
    dict_(dictionary::null),
    parcelTypeId_(parcelTypeId),
    rhoMin_(rhoMin),
    rho0_(rho0),
    minParticleMass_(minParticleMass),
    youngsModulus_(youngsModulus),
    poissonsRatio_(poissonsRatio),
    Cy_(Cy) 
{}


template<class ParcelType>
inline Foam::StreamParcel<ParcelType>::StreamParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    dTarget_(0.0),
    U_(vector::zero),
    f_(vector::zero),
    angularMomentum_(vector::zero),
    torque_(vector::zero),
    rho_(0.0),
    age_(0.0),
    cV_(0.0),
    //- modified by Patrick
    //Mfp_(0.0),   
    //Tp_(300), 
    //YTTIPp_(0.0), 
    //RRp_(0.0),
    N_(1.0),
    A_(0.0),
    V_(0.0),
    a_(0.0),
    v_(0.0),    
    dp_(0.0),
    rc_(0.0),
    cn_(0.0),
    Dp_(0.0),
    Kn_(0.0),
    Cun_(0.0),
    l_(0.0),
    gp_(0.0),
    beta_(0.0),
    tau_(GREAT),        
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0),
    //Mfc_(0.0),
    //epsilonc_(0.0),
    //kc_(0.0),
    //YTTIPc_(0.0)
    Tc_(300.0),
    concSourcec_(0.0),
    gradTc_(vector::zero)
{}


template<class ParcelType>
inline Foam::StreamParcel<ParcelType>::StreamParcel
(
    const polyMesh& owner,
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
    //- modified by Patrick
    //const scalar Mfp0,
    //const scalar YTTIPp0,
    const constantProperties& constProps
)
:
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    active_(true),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    dTarget_(dTarget0),
    U_(U0),
    f_(f0),
    angularMomentum_(angularMomentum0),
    torque_(torque0),
    rho_(constProps.rho0()),
    age_(0.0),
    cV_(0.0),
    //- Modified by Patrick
    //Mfp_(Mfp0),   
    //Tp_(300), 
    //YTTIPp_(YTTIPp0), 
    //RRp_(0.0),
    N_(1.0),
    A_(0.0),
    V_(0.0),
    a_(0.0),
    v_(0.0),    
    dp_(0.0),
    rc_(0.0),
    cn_(0.0),
    Dp_(0.0),
    Kn_(0.0),
    Cun_(0.0),
    l_(0.0),
    gp_(0.0),
    beta_(0.0),
    tau_(GREAT),    
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0),
    //Mfc_(0.0),
    //epsilonc_(0.0),
    //kc_(0.0),
    //YTTIPc_(0.0)    
    //- added by Patrick
    Tc_(300.0),
    concSourcec_(0.0),
    gradTc_(vector::zero)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline const Foam::dictionary&
Foam::StreamParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template<class ParcelType>
inline Foam::label
Foam::StreamParcel<ParcelType>::constantProperties::parcelTypeId() const
{
    return parcelTypeId_;
}


template<class ParcelType>
inline Foam::scalar
Foam::StreamParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_;
}


template<class ParcelType>
inline Foam::scalar
Foam::StreamParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_;
}


template<class ParcelType>
inline Foam::scalar
Foam::StreamParcel<ParcelType>::constantProperties::minParticleMass() const
{
    return minParticleMass_;
}


template<class ParcelType>
inline Foam::scalar
Foam::StreamParcel<ParcelType>::constantProperties::youngsModulus() const
{
    return youngsModulus_;
}


template<class ParcelType>
inline Foam::scalar
Foam::StreamParcel<ParcelType>::constantProperties::poissonsRatio() const
{
    return poissonsRatio_;
}

template<class ParcelType>
inline Foam::scalar
Foam::StreamParcel<ParcelType>::constantProperties::Cy() const
{
    return Cy_;
}

// * * * * * * * StreamParcel Member Functions  * * * * * * * //

template<class ParcelType>
inline bool Foam::StreamParcel<ParcelType>::active() const
{
    return active_;
}


template<class ParcelType>
inline Foam::label Foam::StreamParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::dTarget() const
{
    return dTarget_;
}


template<class ParcelType>
inline const Foam::vector& Foam::StreamParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::StreamParcel<ParcelType>::f() const
{
    return f_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::StreamParcel<ParcelType>::angularMomentum() const
{
    return angularMomentum_;
}


template<class ParcelType>
inline const Foam::vector&
Foam::StreamParcel<ParcelType>::torque() const
{
    return torque_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::age() const
{
    return age_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::cV() const
{
    return cV_;
}

//- modified by Patrick
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::Mfp() const
//{
//    return Mfp_;
//}
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::Tp() const
//{
//    return Tp_;
//}
//
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::YTTIPp() const
//{
//    return YTTIPp_;
//}
//
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::RRp() const
//{
//    return RRp_;
//}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::N() const
{
    return N_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::A() const
{
    return A_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::V() const
{
    return V_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::a() const
{
    return a_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::v() const
{
    return v_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::dp() const
{
    return dp_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::rc() const
{
    return rc_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::cn() const
{
    return cn_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::Dp() const
{
    return Dp_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::Kn() const
{
    return Kn_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::Cun() const
{
    return Cun_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::l() const
{
    return l_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::gp() const
{
    return gp_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::beta() const
{
    return beta_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::tau() const
{
    return tau_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template<class ParcelType>
inline const Foam::vector& Foam::StreamParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::rhoc() const
{
    return rhoc_;
}


template<class ParcelType>
inline const Foam::vector& Foam::StreamParcel<ParcelType>::Uc() const
{
    return Uc_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::muc() const
{
    return muc_;
}

//- modified by Patrick
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::Mfc() const
//{
//    return Mfc_;
//}
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::epsilonc() const
//{
//    return epsilonc_;
//}
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::kc() const
//{
//    return kc_;
//}
//
//
//template<class ParcelType>
//inline Foam::scalar Foam::StreamParcel<ParcelType>::YTTIPc() const
//{
//    return YTTIPc_;
//}

//- added by Patrick

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::Tc() const
{
    return Tc_;
}

template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::concSourcec() const
{
    return concSourcec_;
}

template<class ParcelType>
inline const Foam::vector& Foam::StreamParcel<ParcelType>::gradTc() const
{
    return gradTc_;
}

template<class ParcelType>
inline bool& Foam::StreamParcel<ParcelType>::active()
{
    return active_;
}


template<class ParcelType>
inline Foam::label& Foam::StreamParcel<ParcelType>::typeId()
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::d()
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::dTarget()
{
    return dTarget_;
}


template<class ParcelType>
inline Foam::vector& Foam::StreamParcel<ParcelType>::U()
{
    return U_;
}


template<class ParcelType>
inline Foam::vector& Foam::StreamParcel<ParcelType>::f()
{
    return f_;
}


template<class ParcelType>
inline Foam::vector& Foam::StreamParcel<ParcelType>::angularMomentum()
{
    return angularMomentum_;
}


template<class ParcelType>
inline Foam::vector& Foam::StreamParcel<ParcelType>::torque()
{
    return torque_;
}


template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::rho()
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::age()
{
    return age_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::cV()
{
    return cV_;
}

//- modified by Patrick
//
//template<class ParcelType>
//inline Foam::scalar& Foam::StreamParcel<ParcelType>::Mfp()
//{
//    return Mfp_;
//}
//
//template<class ParcelType>
//inline Foam::scalar& Foam::StreamParcel<ParcelType>::Tp()
//{
//    return Tp_;
//}
//
//template<class ParcelType>
//inline Foam::scalar& Foam::StreamParcel<ParcelType>::YTTIPp()
//{
//    return YTTIPp_;
//}
//
//template<class ParcelType>
//inline Foam::scalar& Foam::StreamParcel<ParcelType>::RRp()
//{
//    return RRp_;
//}


template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::N()
{
    return N_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::A()
{
    return A_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::V()
{
    return V_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::a()
{
    return a_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::v()
{
    return v_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::dp()
{
    return dp_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::rc()
{
    return rc_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::cn()
{
    return cn_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::Dp()
{
    return Dp_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::Kn()
{
    return Kn_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::Cun()
{
    return Cun_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::l()
{
    return l_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::gp()
{
    return gp_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::beta()
{
    return beta_;
}

template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::tau()
{
    return tau_;
}


template<class ParcelType>
inline Foam::scalar& Foam::StreamParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template<class ParcelType>
inline Foam::vector& Foam::StreamParcel<ParcelType>::UTurb()
{
    return UTurb_;
}


template<class ParcelType>
inline Foam::label Foam::StreamParcel<ParcelType>::faceInterpolation() const
{
    // Use volume-based interpolation if dealing with external faces
    if (this->cloud().internalFace(this->face()))
    {
        return this->face();
    }
    else
    {
        return -1;
    }
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::massCell
(
    const label cellI
) const
{
    return rhoc_*this->mesh().cellVolumes()[cellI];
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template<class ParcelType>
inline Foam::vector Foam::StreamParcel<ParcelType>::omega() const
{
    return angularMomentum_/momentOfInertia();
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::areaP(const scalar d)
{
    return 0.25*areaS(d);
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::areaS(const scalar d)
{
    return pi*d*d;
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::Re
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar muc
) const
{
    return rhoc*mag(U - Uc_)*d/(muc + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::We
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar sigma
) const
{
    return rhoc*magSqr(U - Uc_)*d/(sigma + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::StreamParcel<ParcelType>::Eo
(
    const vector& a,
    const scalar d,
    const scalar sigma
) const
{
    vector dir = U_/(mag(U_) + ROOTVSMALL);
    return mag(a & dir)*(rho_ - rhoc_)*sqr(d)/(sigma + ROOTVSMALL);
}


// ************************************************************************* //
