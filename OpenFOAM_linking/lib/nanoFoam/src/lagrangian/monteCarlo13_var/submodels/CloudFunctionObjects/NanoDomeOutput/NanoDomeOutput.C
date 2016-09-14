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

#include "NanoDomeOutput.H"
#include "Pstream.H"
#include "stringListOps.H"
#include "ListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::makeLogFile()
{
        // Create the output file if not already created
        if (debug)
        {
            Info<< "Creating output file" << endl;
        }


   const fileName logDir = outputDir_/("TrackData");

   // Create directory if does not exist

   if (Pstream::master())
   {
     
       mkDir(logDir);

   }
}

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::defineProps()
{
    scalar nPhys = 6;                    // T,P,X,Y,Z = 5

    nSpecs_ = 1;                  // TODO: should read the number of species when real combustion is used
    nProps_ = nPhys + nSpecs_;    // Total number of properties

    propNames_.setSize(nProps_ + nSpecs_);

    propNames_[0] = "TIME_SAMPLES";
    propNames_[1] = "T";
    propNames_[2] = "P";
    propNames_[3] = "X";
    propNames_[4] = "Y";
    propNames_[5] = "Z";
    propNames_[6] = "MOLAR_C";
}

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::extractData(List<string>& sProps, const parcelType& pp)
{
    //List<string> sProps(nProps_);     //Helper Objects for casting
    List<OStringStream> oss(nProps_);
    scalar pressure(1.0e5);           //no pressure on particle at the monemt. TODO: add interpolation etc.

    const scalar simTime = this->owner().time().value();

    //- procedure that directly converts the scalars into strings using the OStringStream of OF
    oss[0] << simTime;	//pp.currentTime();
    sProps[0] = oss[0].str();
    oss[1] << pp.Tc();
    sProps[1] = oss[1].str();
    oss[2] << pressure;
    sProps[2] = oss[2].str();
    oss[3] << pp.position()[0]; 
    sProps[3] = oss[3].str();
    oss[4] << pp.position()[1];
    sProps[4] = oss[4].str();
    oss[5] << pp.position()[2];
    sProps[5] = oss[5].str();
    oss[6] << pp.concSourcec();
    sProps[6] = oss[6].str();
}

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::Create_XML_Node(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement *root_e, string TAG, string text)
{
    using namespace tinyxml2;

    XMLElement *p_new_elem = doc.NewElement(TAG.c_str());
    p_new_elem->SetText(text.c_str());
    root_e->InsertEndChild(p_new_elem);
}

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::writeXMLStream(tinyxml2::XMLDocument& doc, List<string> data, tinyxml2::XMLElement *r_sub_tree, label trackI)
{
	using namespace tinyxml2;

        //- Convert nSpecies_ to a string for XML output
	OStringStream speciesSS;
	speciesSS << nSpecs_;

        //- Convert TrackI to a string for XML output
	OStringStream trackISS;
	trackISS << trackI;

	XMLElement *cursor_ptr = doc.NewElement("STREAM");
	// Write Stream ID
	Create_XML_Node(doc, cursor_ptr, "ID", trackISS.str());
	// Write Time samples
	Create_XML_Node(doc, cursor_ptr, "N_TIME_SAMPLES", nTimes_);
	// Write Time Samples
	Create_XML_Node(doc, cursor_ptr, "TIME_SAMPLES", data[0]);
	// Write Temperatures
	Create_XML_Node(doc, cursor_ptr, "T", data[1]);
	// Write Presssures
	Create_XML_Node(doc, cursor_ptr, "P", data[2]);
	// Write position x
	Create_XML_Node(doc, cursor_ptr, "X", data[3]);
	// Write position y
	Create_XML_Node(doc, cursor_ptr, "Y", data[4]);
	// Write position z
	Create_XML_Node(doc, cursor_ptr, "Z", data[5]);
	// Write number of species defiend as string stream
	Create_XML_Node(doc, cursor_ptr, "N_SPECIES", speciesSS.str());
	// Write Species formula
	Create_XML_Node(doc, cursor_ptr, "SPECIES", "Fe(CO)5");
	//Write Molar Concentrations
	Create_XML_Node(doc, cursor_ptr, "MOLAR_C", data[6]);
	// Add the STREAM Sub-tree to the GP root
	r_sub_tree->InsertEndChild(cursor_ptr);
}


template<class CloudType>
tinyxml2::XMLElement* Foam::NanoDomeOutput<CloudType>::makeXMLTree(tinyxml2::XMLDocument& doc)
{
    using namespace tinyxml2;

    //- Create root node with tag "GP"
    XMLNode *root_pt = doc.NewElement("GP");
    //- Create a moving pointer to element
    XMLElement *cursor_pt = root_pt->ToElement();
    //- Add a new element (child) to the XML tree
    doc.InsertEndChild(cursor_pt);

    //- Create N_STREAMS leaf
    Create_XML_Node(doc, cursor_pt, "N_STREAM", "Number of Strems");

    //- Create T_START leaf 
    Create_XML_Node(doc, cursor_pt, "T_START", "starting sampling time");

    //- Create T_END leaf
    Create_XML_Node(doc, cursor_pt, "T_END", "end sampling time");

    return cursor_pt; 
}

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::saveXML(tinyxml2::XMLDocument& doc)
{
    if(Pstream::master())
    {
        using namespace tinyxml2;

	string logDir;
        
	if(Pstream::parRun())
	{
            logDir = "../Streamline.xml";
	}
	else
	{
            logDir = "Streamline.xml";
	}
        //Save XML FIle in the specified path
        XMLError err = doc.SaveFile(logDir.c_str());
    }
}

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const Time& time = mesh.time();

    const label procI = Pstream::myProcNo();

    Info<< type() << "Writing XML file... ";

    //- list for XML-output listed: tracks< props<string> >       it is a chain of time values
    List<List<string> > XMLData(nTracks_);

    //- TODO: maybe this declaration can be done only for the master node !? 
    //- Define XML document
    tinyxml2::XMLDocument doc;

    //- Create the XML tree and retrun a pointer to the root node
    tinyxml2::XMLElement *cursor_pt = makeXMLTree(doc);

    forAll(particleData_, trackI)
    {
	//- Listed: procs< times<scalar> >
        List<List<scalar> > procTimes(Pstream::nProcs());
	procTimes[Pstream::myProcNo()] = times_[trackI];
        Pstream::gatherList(procTimes);

	//- Listed: procs< times< props<string> > >
	List<List<List<string> > > procData(Pstream::nProcs());
        procData[Pstream::myProcNo()] = particleData_[trackI];
        Pstream::gatherList(procData);

	if(Pstream::master())
	{
	    //- Listed: times< props<string> >    for a single track
	    List<List<string> > globalData;
	    globalData = ListListOps::combine<List<List<string> > >
            (
                procData,
	        accessOp<List<List<string> > >()
            );

            List<scalar> globalTimes;
            globalTimes = ListListOps::combine<List<scalar> >
            (
                procTimes,
                accessOp<List<scalar> >()
            );

	    //- Define nTimes_ for XML output
	    OStringStream tt;
	    tt << globalData.size();
	    nTimes_ = tt.str();

	    //- Find right order of list, which is not in right chronological order
            labelList indices;
            sortedOrder(globalTimes, indices);

	    List<string> tempData(nProps_, "");

	    //- Convert listed data to a chain of time values for NanoDome XML file
	    forAll(globalData, i)
	    {
	        //- Take the right order in time distributed in the order of "indeces"
                label timeI = indices[i];

	        //- Saving the data per time in a new list
	        List<string> timeData;
	        timeData = globalData[timeI];

	        forAll(timeData, propI)
	        {
	    	//Info << propNames_[propI] << nl << timeData[propI] << endl;
                    tempData[propI] += timeData[propI] + " ";
	        }
	    }

	    XMLData[trackI] = tempData;

	    //- Insert data to subtree "streams"
	    writeXMLStream(doc, tempData, cursor_pt, trackI);

	}

    }

    //- Save and close the XML File
    saveXML(doc);

    Info << "done!" << nl;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanoDomeOutput<CloudType>::NanoDomeOutput
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudFunctionObject<CloudType>(dict, owner, typeName),
    parcelType_(this->coeffDict().lookupOrDefault("parcelType", -1)),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    outputFilePtr_(),
    outputDir_(owner.mesh().time().path()),
    times_(),
    particleData_(),
    //- TODO: auto define number of tracks; now only one track, but code can handle more
    nTracks_(4),
    particleId_(),
    procId_(),
    propNames_(),
    nProps_(),
    nSpecs_(),
    nTimes_(-1)
{
    // --- New data for sampling --- //

    // -- TODO : now only for one particle -> multi particle
    particleData_.setSize(nTracks_);
    times_.setSize(nTracks_);
    particleId_.setSize(nTracks_);
    procId_.setSize(nTracks_);
    //-TODO: replace this hardcoded IDs with something chosen on a smart way
    particleId_[0] = 0;
    particleId_[1] = 1; 
    particleId_[2] = 2;
    particleId_[3] = 3;
    procId_[0] = 0;
    procId_[1] = 0;
    procId_[2] = 0;
    procId_[3] = 0;
    //-End of hardcoded stuff
    defineProps();

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

    makeLogFile();

}


template<class CloudType>
Foam::NanoDomeOutput<CloudType>::NanoDomeOutput
(
    const NanoDomeOutput<CloudType>& pc
)
:
    CloudFunctionObject<CloudType>(pc),
    parcelType_(pc.parcelType_),
    resetOnWrite_(pc.resetOnWrite_),
    outputFilePtr_(),
    outputDir_(pc.outputDir_),
    times_(pc.times_),
    particleData_(pc.particleData_),
    nTracks_(pc.nTracks_),
    particleId_(pc.particleId_),
    procId_(pc.procId_),
    propNames_(pc.propNames_),
    nProps_(pc.nProps_),
    nSpecs_(pc.nSpecs_),
    nTimes_(pc.nTimes_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NanoDomeOutput<CloudType>::~NanoDomeOutput()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NanoDomeOutput<CloudType>::nanoDomeCall
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


    forAll(particleData_, trackI)  //Maybe use a map instead of an loop //  TODO: insert map
    {
       if (labelPair(p.origProc(), p.origId()) == labelPair(procId_[trackI], particleId_[trackI]))
       {

           times_[trackI].append(this->owner().time().value());

	   List<string> returnData(nProps_);
           extractData(returnData, p);
	   particleData_[trackI].append(returnData);

       }
    }
}


// ************************************************************************* //
