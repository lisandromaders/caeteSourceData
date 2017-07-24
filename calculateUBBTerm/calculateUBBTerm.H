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


\*---------------------------------------------------------------------------*/

#ifndef calculateUBBTerm_H
#define calculateUBBTerm_H

//#include "basicPsiThermo.H"
#include "readUBBData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class calculateUBBTerm Declaration
\*---------------------------------------------------------------------------*/

class calculateUBBTerm
{
	    //- Table look up routines
        const readUBBData& UBBTable_;
        
                
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        calculateUBBTerm(const calculateUBBTerm&);


public:

    //- Runtime type information
    TypeName("calculateUBBTerm");


    // Constructors

        //- Construct from mesh
        calculateUBBTerm
        (
            const fvMesh&,
            const readUBBData&
			//sigmaMHD
        );


    //- Destructor
    ~calculateUBBTerm();


    // Member functions

        //- Update properties
        void correct();

	virtual tmp<volVectorField> UBB() const
	{
		return UBB_;
	}

	virtual volVectorField& UBB()
	{
		return UBB_;
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// #include "hPsiFanzyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
