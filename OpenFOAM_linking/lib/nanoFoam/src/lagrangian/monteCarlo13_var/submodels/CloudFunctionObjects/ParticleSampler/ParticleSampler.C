/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "ParticleSampler.H"
#include "Pstream.H"
#include "stringListOps.H"
#include "ListOps.H"
#include "ListListOps.H"

#include "surfaceWriter.H"
#include "unitConversion.H"
#include "Random.H"
#include "triangle.H"
#include "cloud.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleSampler<CloudType>::makeLogFile
(
    const faceList& faces,
    const Field<point>& points,
    const Field<scalar>& area
)
{
    // Create the output file if not already created
    if (log_)
    {
        if (debug)
        {
            Info<< "Creating output file" << endl;
        }

        if (Pstream::master())
        {
            const fileName logDir = outputDir_/this->owner().time().timeName();

            // Create directory if does not exist
            mkDir(logDir);

            // Open new file at start up
            outputFilePtr_.reset
            (
                new OFstream(logDir/(type() + ".dat"))
            );

            outputFilePtr_()
                << "# Source     : " << type() << nl
                << "# Total area : " << sum(area) << nl;

            //forAll(faces, i)
            //{
            //    word id = Foam::name(i);

            //    outputFilePtr_()
            //        << tab << "area[" << id << "]"
            //        << tab << "mass[" << id << "]"
            //        << tab << "massFlowRate[" << id << "]"
            //        << endl;
            //}
	   
	    // ---- Writing header -.- New ---- //  
            string header("# Time currentProc " + parcelType::propertyList_);
            outputFilePtr_() << header.c_str() << nl;

        }
    }
}


template<class CloudType>
void Foam::ParticleSampler<CloudType>::initPolygons()
{
    mode_ = mtPolygon;

    List<Field<point> > polygons(this->coeffDict().lookup("polygons"));
    label nPoints = 0;
    forAll(polygons, polyI)
    {
        label np = polygons[polyI].size();
        if (np < 3)
        {
            FatalIOErrorIn
            (
                "Foam::ParticleSampler<CloudType>::initPolygons()",
                this->coeffDict()
            )
                << "polygons must consist of at least 3 points"
                << exit(FatalIOError);
        }

        nPoints += np;
    }

    label pointOffset = 0;
    points_.setSize(nPoints);
    faces_.setSize(polygons.size());
    faceTris_.setSize(polygons.size());
    area_.setSize(polygons.size());
    forAll(faces_, faceI)
    {
        const Field<point>& polyPoints = polygons[faceI];
        face f(identity(polyPoints.size()) + pointOffset);
        UIndirectList<point>(points_, f) = polyPoints;
        area_[faceI] = f.mag(points_);

        DynamicList<face> tris;
        f.triangles(points_, tris);
        faceTris_[faceI].transfer(tris);

        faces_[faceI].transfer(f);

        pointOffset += polyPoints.size();
    }
}


template<class CloudType>
void Foam::ParticleSampler<CloudType>::initConcentricCircles()
{
    mode_ = mtConcentricCircle;

    vector origin(this->coeffDict().lookup("origin"));

    radius_ = this->coeffDict().lookup("radius");
    nSector_ = readLabel(this->coeffDict().lookup("nSector"));

    label nS = nSector_;

    vector refDir;
    if (nSector_ > 1)
    {
        refDir = this->coeffDict().lookup("refDir");
        refDir -= normal_*(normal_ & refDir);
        refDir /= mag(refDir);
    }
    else
    {
        // set 4 quadrants for single sector cases
        nS = 4;

        vector tangent = vector::zero;
        scalar magTangent = 0.0;

        Random rnd(1234);
        while (magTangent < SMALL)
        {
            vector v = rnd.vector01();

            tangent = v - (v & normal_)*normal_;
            magTangent = mag(tangent);
        }

        refDir = tangent/magTangent;
    }

    scalar dTheta = 5.0;
    scalar dThetaSector = 360.0/scalar(nS);
    label intervalPerSector = max(1, ceil(dThetaSector/dTheta));
    dTheta = dThetaSector/scalar(intervalPerSector);

    label nPointPerSector = intervalPerSector + 1;

    label nPointPerRadius = nS*(nPointPerSector - 1);
    label nPoint = radius_.size()*nPointPerRadius;
    label nFace = radius_.size()*nS;

    // add origin
    nPoint++;

    points_.setSize(nPoint);
    faces_.setSize(nFace);
    area_.setSize(nFace);

    coordSys_ = cylindricalCS("coordSys", origin, normal_, refDir, false);

    List<label> ptIDs(identity(nPointPerRadius));

    points_[0] = origin;

    // points
    forAll(radius_, radI)
    {
        label pointOffset = radI*nPointPerRadius + 1;

        for (label i = 0; i < nPointPerRadius; i++)
        {
            label pI = i + pointOffset;
            point pCyl(radius_[radI], degToRad(i*dTheta), 0.0);
            points_[pI] = coordSys_.globalPosition(pCyl);
        }
    }

    // faces
    DynamicList<label> facePts(2*nPointPerSector);
    forAll(radius_, radI)
    {
        if (radI == 0)
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                // append origin point
                facePts.append(0);

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = ptIDs.fcIndex(i - 1) + 1;
                    facePts.append(id);
                }

                label faceI = secI + radI*nS;

                faces_[faceI] = face(facePts);
                area_[faceI] = faces_[faceI].mag(points_);
            }
        }
        else
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                label offset = (radI - 1)*nPointPerRadius + 1;

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }
                for (label ptI = nPointPerSector-1; ptI >= 0; ptI--)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + nPointPerRadius + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }

                label faceI = secI + radI*nS;

                faces_[faceI] = face(facePts);
                area_[faceI] = faces_[faceI].mag(points_);
            }
        }
    }
}


template<class CloudType>
Foam::label Foam::ParticleSampler<CloudType>::collectParcelPolygon
(
    const point& position,
    const vector& U
) const
{
    scalar dt = this->owner().db().time().deltaTValue();

    point end(position + dt*U);

    label dummyNearType = -1;
    label dummyNearLabel = -1;

    forAll(faces_, faceI)
    {
        const label facePoint0 = faces_[faceI][0];

        const point p0 = points_[facePoint0];

        const scalar d1 = normal_ & (position - p0);
        const scalar d2 = normal_ & (end - p0);

        if (sign(d1) == sign(d2))
        {
            // did not cross polygon plane
            continue;
        }

        // intersection point
        point pIntersect = position + (d1/(d1 - d2))*dt*U;

        const List<face>& tris = faceTris_[faceI];

        // identify if point is within poly bounds
        forAll(tris, triI)
        {
            const face& tri = tris[triI];
            triPointRef t
            (
                points_[tri[0]],
                points_[tri[1]],
                points_[tri[2]]
            );

            if (t.classify(pIntersect, dummyNearType, dummyNearLabel))
            {
                return faceI;
            }
        }
    }

    return -1;
}


template<class CloudType>
Foam::label Foam::ParticleSampler<CloudType>::collectParcelConcentricCircles
(
    const point& position,
    const vector& U
) const
{
    label secI = -1;

    scalar dt = this->owner().db().time().deltaTValue();

    point end(position + dt*U);

    const scalar d1 = normal_ & (position - coordSys_.origin());
    const scalar d2 = normal_ & (end - coordSys_.origin());

    if (sign(d1) == sign(d2))
    {
        // did not cross plane
        return secI;
    }

    // intersection point in cylindrical co-ordinate system
    point pCyl = coordSys_.localPosition(position + (d1/(d1 - d2))*dt*U);

    scalar r = pCyl[0];

    if (r < radius_.last())
    {
        label radI = 0;
        while (r > radius_[radI])
        {
            radI++;
        }

        if (nSector_ == 1)
        {
            secI = 4*radI;
        }
        else
        {
            scalar theta = pCyl[1] + constant::mathematical::pi;

            secI =
                nSector_*radI
              + floor
                (
                    scalar(nSector_)*theta/constant::mathematical::twoPi
                );
        }
    }

    return secI;
}

//template<class CloudType>
//void Foam::ParticleSampler<CloudType>::write()
//{
//    forAll(faceData_, i)
//    {
//        List<List<scalar> > procTimes(Pstream::nProcs());
//        procTimes[Pstream::myProcNo()] = times_[i];
//        Pstream::gatherList(procTimes);
//
//        List<List<string> > procData(Pstream::nProcs());
//        procData[Pstream::myProcNo()] = faceData_[i];
//        Pstream::gatherList(procData);
//
//        if (Pstream::master())
//        {
//            const fvMesh& mesh = this->owner().mesh();
//
//            fileName outputDir = mesh.time().path();
//
//            if (Pstream::parRun())
//            {
//                // Put in undecomposed case (Note: gives problems for
//                // distributed data running)
//                outputDir =
//                    outputDir/".."/"postProcessing"/cloud::prefix/
//                    this->owner().name()/mesh.time().timeName();
//            }
//            else
//            {
//                outputDir =
//                    outputDir/"postProcessing"/cloud::prefix/
//                    this->owner().name()/mesh.time().timeName();
//            }
//
//            // Create directory if it doesn't exist
//            mkDir(outputDir);
//
//	    word id = Foam::name(i);
//
//            const word& faceName = "Poly_"+id;
//
//            OFstream faceOutFile
//            (
//                outputDir/faceName + ".post",
//                IOstream::ASCII,
//                IOstream::currentVersion,
//                mesh.time().writeCompression()
//            );
//
//            List<string> globalData;
//            globalData = ListListOps::combine<List<string> >
//            (
//                procData,
//                accessOp<List<string> >()
//            );
//
//            List<scalar> globalTimes;
//            globalTimes = ListListOps::combine<List<scalar> >
//            (
//                procTimes,
//                accessOp<List<scalar> >()
//            );
//
//            labelList indices;
//            sortedOrder(globalTimes, indices);
//
//            string header("# Time currentProc " + parcelType::propertyList_);
//            faceOutFile<< header.c_str() << nl;
//
//            forAll(globalTimes, i)
//            {
//                label dataI = indices[i];
//
//                faceOutFile
//                    << globalTimes[dataI] << ' '
//                    << globalData[dataI].c_str()
//                    << nl;
//            }
//        }
//
//        faceData_[i].clearStorage();
//        times_[i].clearStorage();
//    }
//}


template<class CloudType>
void Foam::ParticleSampler<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const Time& time = mesh.time();

    const label procI = Pstream::myProcNo();

    Info<< type() << " output:" << nl;

    //if (outputFilePtr_.valid())
    //{
    //    outputFilePtr_() << time.timeName();
    //}

    //Field<scalar> faceMassTotal(mass_.size(), 0.0);
    //this->getModelProperty("massTotal", faceMassTotal);

    //Field<scalar> faceMassFlowRate(massFlowRate_.size(), 0.0);
    //this->getModelProperty("massFlowRate", faceMassFlowRate);

    forAll(faces_, faceI)
    {
        List<List<scalar> > procTimes(Pstream::nProcs());
        procTimes[Pstream::myProcNo()] = times_[faceI];
        Pstream::gatherList(procTimes);

        List<List<string> > procData(Pstream::nProcs());
        procData[Pstream::myProcNo()] = faceData_[faceI];
        Pstream::gatherList(procData);

        List<string> globalData;
        globalData = ListListOps::combine<List<string> >
        (
            procData,
            accessOp<List<string> >()
        );

        List<scalar> globalTimes;
        globalTimes = ListListOps::combine<List<scalar> >
        (
            procTimes,
            accessOp<List<scalar> >()
        );


        labelList indices;
        sortedOrder(globalTimes, indices);

        if (outputFilePtr_.valid())
        {
            forAll(globalTimes, i)
            {
                label dataI = indices[i];

                outputFilePtr_()
                    << globalTimes[dataI] << ' '
                    << globalData[dataI].c_str()
                    << nl;
            }

            //outputFilePtr_()
            //    << tab << area_[faceI]
            //    << tab << faceMassTotal[faceI]
            //    << tab << faceMassFlowRate[faceI]
            //    << endl;
        }
        faceData_[faceI].clearStorage();
        times_[faceI].clearStorage();
    }

    //if (surfaceFormat_ != "none")
    //{
    //    if (Pstream::master())
    //    {
    //        autoPtr<surfaceWriter> writer(surfaceWriter::New(surfaceFormat_));

    //        writer->write
    //        (
    //            outputDir_/time.timeName(),
    //            "collector",
    //            points_,
    //            faces_,
    //            "massTotal",
    //            faceMassTotal,
    //            false
    //        );

    //        writer->write
    //        (
    //            outputDir_/time.timeName(),
    //            "collector",
    //            points_,
    //            faces_,
    //            "massFlowRate",
    //            faceMassFlowRate,
    //            false
    //        );
    //    }
    //}


    //if (resetOnWrite_)
    //{
    //    Field<scalar> dummy(faceMassTotal.size(), 0.0);
    //    this->setModelProperty("massTotal", dummy);
    //    this->setModelProperty("massFlowRate", dummy);

    //    timeOld_ = timeNew;
    //    totalTime_ = 0.0;
    //}
    //else
    //{
    //    this->setModelProperty("massTotal", faceMassTotal);
    //    this->setModelProperty("massFlowRate", faceMassFlowRate);
    //}

    //forAll(faces_, faceI)
    //{
    //    mass_[faceI] = 0.0;
    //    massTotal_[faceI] = 0.0;
    //    massFlowRate_[faceI] = 0.0;
    //}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleSampler<CloudType>::ParticleSampler
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudFunctionObject<CloudType>(dict, owner, typeName),
    mode_(mtUnknown),
    parcelType_(this->coeffDict().lookupOrDefault("parcelType", -1)),
    removeCollected_(this->coeffDict().lookup("removeCollected")),
    points_(),
    faces_(),
    faceTris_(),
    nSector_(0),
    radius_(),
    coordSys_(false),
    normal_(this->coeffDict().lookup("normal")),
    negateParcelsOppositeNormal_
    (
        readBool(this->coeffDict().lookup("negateParcelsOppositeNormal"))
    ),
    surfaceFormat_(this->coeffDict().lookup("surfaceFormat")),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    totalTime_(0.0),
    mass_(),
    massTotal_(),
    massFlowRate_(),
    log_(this->coeffDict().lookup("log")),
    outputFilePtr_(),
    outputDir_(owner.mesh().time().path()),
    timeOld_(owner.mesh().time().value()),
    times_(),
    faceData_()
{
    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        outputDir_ =
            outputDir_/".."/"postProcessing"/cloud::prefix/owner.name();
    }
    else
    {
        outputDir_ =
            outputDir_/"postProcessing"/cloud::prefix/owner.name();
    }

    normal_ /= mag(normal_);

    word mode(this->coeffDict().lookup("mode"));
    if (mode == "polygon")
    {
        initPolygons();
    }
    else if (mode == "concentricCircle")
    {
        initConcentricCircles();
    }
    else
    {
        FatalErrorIn
        (
            "Foam::ParticleSampler<CloudType>::ParticleSampler"
            "("
                "const dictionary&,"
                "CloudType&"
            ")"
        )
            << "Unknown mode " << mode << ".  Available options are "
            << "polygon and concentricCircle" << exit(FatalError);
    }

    mass_.setSize(faces_.size(), 0.0);
    massTotal_.setSize(faces_.size(), 0.0);
    massFlowRate_.setSize(faces_.size(), 0.0);

    makeLogFile(faces_, points_, area_);

    // --- New data for sampling --- //

    faceData_.setSize(faces_.size());
    times_.setSize(faces_.size());
}


template<class CloudType>
Foam::ParticleSampler<CloudType>::ParticleSampler
(
    const ParticleSampler<CloudType>& pc
)
:
    CloudFunctionObject<CloudType>(pc),
    mode_(pc.mode_),
    parcelType_(pc.parcelType_),
    removeCollected_(pc.removeCollected_),
    points_(pc.points_),
    faces_(pc.faces_),
    faceTris_(pc.faceTris_),
    nSector_(pc.nSector_),
    radius_(pc.radius_),
    coordSys_(pc.coordSys_),
    normal_(pc.normal_),
    negateParcelsOppositeNormal_(pc.negateParcelsOppositeNormal_),
    surfaceFormat_(pc.surfaceFormat_),
    resetOnWrite_(pc.resetOnWrite_),
    totalTime_(pc.totalTime_),
    mass_(pc.mass_),
    massTotal_(pc.massTotal_),
    massFlowRate_(pc.massFlowRate_),
    log_(pc.log_),
    outputFilePtr_(),
    outputDir_(pc.outputDir_),
    timeOld_(0.0),
    times_(pc.times_),
    faceData_(pc.faceData_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleSampler<CloudType>::~ParticleSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleSampler<CloudType>::postMove
(
    const parcelType& p,
    const label cellI,
    const scalar dt,
    bool& keepParticle
)
{
    if ((parcelType_ != -1) && (parcelType_ != p.typeId()))
    {
        return;
    }

    label faceI = -1;

    switch (mode_)
    {
        case mtPolygon:
        {
            faceI = collectParcelPolygon(p.position(), p.U());
            break;
        }
        case mtConcentricCircle:
        {
            faceI = collectParcelConcentricCircles(p.position(), p.U());
            break;
        }
        default:
        {
        }
    }

    if (faceI != -1)
    {

	times_[faceI].append(this->owner().time().value());

	OStringStream data;
	data<< Pstream::myProcNo() << ' ' << p;

	faceData_[faceI].append(data.str());

	// --------- Outcommented the mass calculation since nanoparticle mass is not properly considered
	
        //scalar m = p.nParticle()*p.mass();

        //if (negateParcelsOppositeNormal_)
        //{
        //    vector Uhat = p.U();
        //    Uhat /= mag(Uhat) + ROOTVSMALL;
        //    if ((Uhat & normal_) < 0)
        //    {
        //        m *= -1.0;
        //    }
        //}

        //// add mass contribution
        //mass_[faceI] += m;

        //if (nSector_ == 1)
        //{
        //    mass_[faceI + 1] += m;
        //    mass_[faceI + 2] += m;
        //    mass_[faceI + 3] += m;
        //}

        if (removeCollected_)
        {
            keepParticle = false;
        }
    }
}


// ************************************************************************* //
