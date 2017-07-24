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

#include "calculateUBBTerm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calculateUBBTerm, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::calculateUBBTerm::calculate(fvVectorMatrix& UEqn)
{
	vectorField& BCells 		= this->B_.internalField();
    vectorField& UBBCells 		= this->UBB_.internalField();
	vectorField& UCells 		= this->UEqn.psi().internalField();
    //scalarField& sigmaMHDCells	= this->sigmaMHD_.internalField();
    
    forAll(UCells, celli)
    {
		BCells[celli] = UBBTable_.assignUBBValues(mesh);
		//UBBCells[celli] = sigmaMHDCells[celli]*(UCells[celli] ^ BCells[celli]) ^ BCells[celli];
		UBBCells[celli] = (UCells[celli] ^ BCells[celli]) ^ BCells[celli];
    }
/* 
    forAll(U_.boundaryField(), patchi)
    {
		fvPatchVectorField& pB 			= this->B_.boundaryField()[patchi];
		fvPatchVectorField& pUBB 		= this->UBB_.boundaryField()[patchi];
		//fvPatchVectorField& pU 			= this->U_.boundaryField()[patchi];
		//fvPatchScalarField& psigmaMHD	= this->sigmaMHD_.boundaryField()[patchi];
        
        if (pB.fixesValue())
        {
            forAll(pUBB, facei)
            {
                pB[facei] = UBBTable_.assignUBBValues(mesh);
				//UBBCells[facei] = sigmaMHDCells[facei]*(UCells[facei] ^ BCells[facei]) ^ BCells[facei];
				UBBCells[facei] = (UCells[facei] ^ BCells[facei]) ^ BCells[facei];
            }
        }
        else
        {
            forAll(pUBB, facei)
            {
                pB[facei] = UBBTable_.assignUBBValues(mesh);
				//UBBCells[facei] = sigmaMHDCells[facei]*(UCells[facei] ^ BCells[facei]) ^ BCells[facei];
				UBBCells[facei] = (UCells[facei] ^ BCells[facei]) ^ BCells[facei];                
            }
        }
    } */
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calculateUBBTerm::calculateUBBTerm
(
    const fvMesh& mesh,
    const readUBBData& readUBBData
    //const readJBData& readJBData,
    //const volScalarField& sigmaMHD
)

    volVectorField B_
    (
        IOobject
        (
            "B",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
		dimensionSet(1,0,-2,0,0,-1,0),
		readUBBData.assignUBBValues(mesh)
    ),
	
	volVectorField UBB_
    (
        IOobject
        (
            "UBB",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
		dimensionSet(2,1,-5,0,0,-2,0)
    ),

    JBTable_(readJBData),
	UBBTable_(readUBBData)
    //sigmaMHD_(sigmaMHD),
    
/* 	//AQUI EU INICIALIZO MEUS VALORES NO PRÓPRIO CONSTRUTOR (o que eu estava fazendo antes no createSourceTerm.H)
{
    vectorField& UBBCells = this->UBB_.internalField();
	//vectorField& UCells = this->U_.internalField(); // acho que eu nao preciso declarar U, ele deve estar em "mesh"

    forAll(UBBCells, celli)
    {
        BCells[celli] = UBBTable_.assignUBBValues(mesh);
		UBBCells[celli] = (UCells[celli] ^ BCells[celli]) ^ BCells[celli];
    }

    forAll(UBB_.boundaryField(), patchi)
    {
        fvPatchVectorField& pUBB 	= this->UBB_.boundaryField()[patchi];
		fvPatchVectorField& pB 		= this->B_.boundaryField()[patchi];
        //fvPatchVectorField& pU 		= this->U_.boundaryField()[patchi];
		
        // Initialize enthalpy and temperature with values from the FGM table
        if (pUBB.fixesValue())
        {
            forAll(pUBB, facei)
            {
                pUBB[facei] = (pU[facei] ^ pB[facei]) ^ pB[facei];
            }
        }
        else
        {
            forAll(pUBB, facei)
            {
                pUBB[facei] = (pU[facei] ^ pB[facei]) ^ pB[facei];
            }
        }
    }

    // Correct enthalpy and temperature boundary conditions after initilisation
    UBB_.correctBoundaryConditions();

/*   
    // With the changes made in the fanzyLookUp class the following lines did not work anymore. Probably, the function dataNames() was created
    // to deal with the old array dimensions.. When everything is compiing and running fine I will take a look on this.
 
   Info << "dataFieldName[rho] = "
        << fgmTable_.dataNames()[fgmThermoTransportIndices_[0]] << nl
        << "dataFieldName[T] = "
        << fgmTable_.dataNames()[fgmThermoTransportIndices_[1]] << nl
        << "dataFieldName[Cp] = "
        << fgmTable_.dataNames()[fgmThermoTransportIndices_[2]] << nl
        << "dataFieldName[kappa] = "
        << fgmTable_.dataNames()[fgmThermoTransportIndices_[3]] << nl
        << "dataFieldName[mu] = "
        << fgmTable_.dataNames()[fgmThermoTransportIndices_[4]] << nl
        << endl;
*/    
    calculate();

    // Switch on saving old time
/*     this->psi_.oldTime();
    this->rho_.oldTime();
    this->mu_.oldTime(); */
	
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calculateUBBTerm::~calculateUBBTerm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calculateUBBTerm::correct()
{
    if (debug)
    {
        Info<< "entering calculateUBBTerm::correct()" << endl;
    }

    // force the saving of the old-time values
/*     this->psi_.oldTime();
    this->rho_.oldTime();
    this->mu_.oldTime(); */

    calculate();

    if (debug)
    {
        Info<< "exiting hPsiFanzy::correct()" << endl;
    }
}

// ************************************************************************* //
