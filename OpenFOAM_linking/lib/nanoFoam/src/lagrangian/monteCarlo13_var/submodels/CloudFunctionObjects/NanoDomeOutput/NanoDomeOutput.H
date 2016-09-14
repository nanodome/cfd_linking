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
    Foam::NanoDomeOutput

Description
    Function object to collect the parcel mass- and mass flow rate over a
    set of polygons.  The polygons are defined as lists of points.  If a
    parcel is 'collected', it is subsequently flagged to be removed from the
    domain.

SourceFiles
    NanoDomeOutput.C

\*---------------------------------------------------------------------------*/

#ifndef NanoDomeOutput_H
#define NanoDomeOutput_H

#include "CloudFunctionObject.H"
#include "Switch.H"
#include "OFstream.H"
#include "tinyxml2.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NanoDomeOutput Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NanoDomeOutput
:
    public CloudFunctionObject<CloudType>
{

private:

    // Private Data

        // Typedefs

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Index of parcel types to collect (-1 by default = all particles)
        const label parcelType_;

        //- Flag to indicate whether data should be reset/cleared on writing
        Switch resetOnWrite_;

        //- Output file pointer
        //autoPtr<OFstream> outputFilePtr_;
	autoPtr<OFstream> outputFilePtr_; 

        //- Output directory
        fileName outputDir_;

	// ----------- New Part ------------- //
	
	//- List of time for each data record listed: tracks< times<scalar> >
	List<DynamicList<scalar> > times_;

	//- List of output data of particles listed: tracks< times< props<string> > >
	//List<DynamicList<string> > particleData_;
	List<DynamicList<List<string> > > particleData_;

	//- Number of particles to track
	const label nTracks_;

	//- orgId of particle which are tracked
	List<label> particleId_;

	//- orgProcId of particles which are tracked
	List<label> procId_;

	//- Names of properties used by XML file
	List<string> propNames_;

	//- Number of properties which are written out
	scalar nProps_;

	//- Number of species which are written out
	scalar nSpecs_;

	//- Numer of time samples
	string nTimes_;

    // Private Member Functions

        //- Helper function to create log files
        void makeLogFile();

	//- Function to create the amound of properties
	void defineProps();

	//- Function to create the amound of properties
	void extractData(List<string>& sProps, const parcelType& pp);

	//- Function to save the XML file
	void saveXML(tinyxml2::XMLDocument& doc);

	//- Function to create the XML tree
	tinyxml2::XMLElement* makeXMLTree(tinyxml2::XMLDocument& doc);

	//- Function to create the XML node
	void Create_XML_Node(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement *root_e, string TAG, string text);

	//- Function to write streams to XML file
	void writeXMLStream(tinyxml2::XMLDocument& doc, List<string> data, tinyxml2::XMLElement *r_sub_tree, label trackI);

protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("nanoDomeOutput");


    // Constructors

        //- Construct from dictionary
        NanoDomeOutput(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NanoDomeOutput(const NanoDomeOutput<CloudType>& pc);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new NanoDomeOutput<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NanoDomeOutput();


    // Member Functions

        // Access

            //- Return const access to the reset on write flag
            inline const Switch& resetOnWrite() const;


        // Evaluation

            //- Post-move hook
            virtual void nanoDomeCall
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

#include "NanoDomeOutputI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NanoDomeOutput.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
