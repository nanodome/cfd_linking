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

template<class ParcelType>
template<class CloudType>
inline Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part
)
:
    ParcelType::template TrackingData<CloudType>(cloud),
    rhoInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.rho()
        )
    ),
    UInterp_
    (
        interpolation<vector>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.U()
        )
    ),
    muInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.mu()
        )
    ),
    //- modified by Patrick
    TInterp_
    (
     	interpolation<scalar>::New
	(
	    cloud.solution().interpolationSchemes(),
	    cloud.T()
	)
    ),
    concSourceInterp_
    (
     	interpolation<scalar>::New
	(
	    cloud.solution().interpolationSchemes(),
	    cloud.concSource()
	)
    ),
    gradTInterp_
    (
     	interpolation<vector>::New
	(
	    cloud.solution().interpolationSchemes(),
	    cloud.gradT()
	)
    ),
    //MfInterp_
    //(
    //    interpolation<scalar>::New
    //    (
    //        cloud.solution().interpolationSchemes(),
    //        cloud.Mf()
    //    )
    //),
    //epsilonInterp_
    //(
    //    interpolation<scalar>::New
    //    (
    //        cloud.solution().interpolationSchemes(),
    //        cloud.epsilon()
    //    )
    //),
    //
    //kInterp_
    //(
    //    interpolation<scalar>::New
    //    (
    //        cloud.solution().interpolationSchemes(),
    //        cloud.k()
    //    )
    //), 
    //
    //YTTIPInterp_
    //(
    //    interpolation<scalar>::New
    //    (
    //        cloud.solution().interpolationSchemes(),
    //        cloud.YTTIP()
    //    )
    //), 
    g_(cloud.g().value()),
    part_(part)
{}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::rhoInterp() const
{
    return rhoInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::UInterp() const
{
    return UInterp_();
}


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::muInterp() const
{
    return muInterp_();
}

//-modified by Patrick 

template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::TInterp() const
{
    return TInterp_();
}

template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::concSourceInterp() const
{
    return concSourceInterp_();
}

template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::gradTInterp() const
{
    return gradTInterp_();
}

//template<class ParcelType>
//template<class CloudType>
//inline const Foam::interpolation<Foam::scalar>&
//Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::MfInterp() const
//{
//    return MfInterp_();
//}
//
//template<class ParcelType>
//template<class CloudType>
//inline const Foam::interpolation<Foam::scalar>&
//Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::epsilonInterp() const
//{
//    return epsilonInterp_();
//}
//
//template<class ParcelType>
//template<class CloudType>
//inline const Foam::interpolation<Foam::scalar>&
//Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::kInterp() const
//{
//    return kInterp_();
//}
//
//template<class ParcelType>
//template<class CloudType>
//inline const Foam::interpolation<Foam::scalar>&
//Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::YTTIPInterp() const
//{
//    return YTTIPInterp_();
//}

template<class ParcelType>
template<class CloudType>
inline const Foam::vector&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::g() const
{
    return g_;
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::StreamParcel<ParcelType>::template
TrackingData<CloudType>::trackPart
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::part() const
{
    return part_;
}


template<class ParcelType>
template<class CloudType>
inline typename Foam::StreamParcel<ParcelType>::template
TrackingData<CloudType>::trackPart&
Foam::StreamParcel<ParcelType>::TrackingData<CloudType>::part()
{
    return part_;
}


// ************************************************************************* //
