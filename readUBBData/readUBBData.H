/*---------------------------------------------------------------------------*\
  
 ####################                                   ###################
 #                  #                                   #                 #   
 #       /||=====   #                                   #   /||  |||||||  #
 #      //||        #                                   #  //||       ||  #
 #     // ||        #     The FanzyFgm tool             #    ||       ||  #
 #    //==||===     #                                   #    ||    |||||  #
 #   //   ||        #     Copyright (C) 2013 by a.f.    #    ||       ||  #
 #  //    ||anzy    #                                   #    ||       ||  #
 #                  #                                   #    ||   ||||||  #
 ####################                                   ===================

-------------------------------------------------------------------------------
License
    This file is part of the CaeteMHD library.

       
Class
    Foam::readUBBData

Description
    Class for the source data looking-up

SourceFiles
    readUBBData.C

\*---------------------------------------------------------------------------*/

#ifndef readUBBData_H
#define readUBBData_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalar.H"
#include "scalarField.H"
#include "volFields.H"
//#include "PtrList.T.H"
#include "PtrList.H"
#include "zeroGradientFvPatchFields.H"
#include "OFstream.H"
#include "IFstream.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class readUBBData Declaration
\*---------------------------------------------------------------------------*/

class readUBBData
:
    public IOdictionary
{
    // Private Member Functions

        //- Disallow copy constructor
        readUBBData(const readUBBData&);

        //- Disallow default bitwise assignment
        void operator=(const readUBBData&);

protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Filename of the FGM table
        const fileName sourceFile_;

      // Variables read from the FGM table
        
          label nPoints_;
          label nColunas_;

          //- FGM data table
          List<scalar> sizeColunas_;

          //- FGM data table
          List<List<scalar> > temp_;

          //- FGM data table
          List<List<scalar> >  sourceData_;

    // Protected Member Functions

        //- Read and discard to newline
        static void readToNewline(IFstream&);

        //- Read FGM table from file
        void readSourceFile(const fileName);

public:

    //- Runtime type information
    TypeName("readUBBData");


    // Constructors

        //- Construct from components
        readUBBData(const fvMesh&);


    //- Destructor
    ~readUBBData();


    // Access general source table data

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        tmp<volVectorField> assignUBBValues(const fvMesh&);
//        void assignSourceValues(const fvMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
