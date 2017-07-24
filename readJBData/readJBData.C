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

       
\*---------------------------------------------------------------------------*/

#include "readJBData.H"
#include "fvMesh.H"
//#include "Time.T.H"
//#include "ListOps.T.H"
#include "Time.H"
#include "ListOps.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(readJBData, 0);
}

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * * //

Foam::readJBData::readJBData
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "mhdProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
	
    mesh_(mesh),
    sourceFile_(fileName(lookup("JB_file")).expand())
{
    Info<< "Iniciando construtor readJBData: arquivo " << sourceFile_ << endl;
    
    readSourceFile(sourceFile_);
}


// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * * //


Foam::readJBData::~readJBData()
{}


// * * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * * * * //

void Foam::readJBData::readToNewline(IFstream& is)
{
    char ch = '\n';
    do
    {
        (is).get(ch);
    }
    while ((is) && ch != '\n');
}

void Foam::readJBData::readSourceFile(const fileName inputName)
{
    IFstream is(inputName);
    //char c;

    if (is.good())
    {
        Info<< "readElectricField: reading forces field file..." << endl;
    }
    else
    {
        FatalErrorIn
        (
            "nomeClasse::readMHDfile(const filename inputName)"
        )
        << "cannot read " << is.name() << abort(FatalError);
    }
    
    nColunas_ = 7;
    sizeColunas_ = List<scalar>(nColunas_);

    // --------  Vamos armazenar os dados numa matriz
    
    temp_ = List<List<scalar> > 
                (
                    1,
                    List<scalar> 
                    (
                        sizeColunas_
                    )
                );
    
    forAll(sizeColunas_,i)
    {
        is.read(temp_[0][i]);
    }
    
    nPoints_ = temp_[0][0];
    
    sourceData_ = List<List<scalar> >
    		(
    		    nPoints_,
    		    List<scalar>
    		    (
    			 sizeColunas_
    		    )
    		);
    
    forAll(sizeColunas_,i)
    {
        sourceData_[0][i] = temp_[0][i];
    }
    
    for (label i=1; i<nPoints_; i++)
    {
        for (label j=0; j<7; j++)
        {
    	is.read(sourceData_[i][j]);
        }
    }
    
    // Mandando imprimir na tela para testar se esta lendo tudo certinho..
    
//   Info<< "Verificando as primeiras 5 linhas para ver se o arquivo esta sendo lido corretamente" << endl << endl;
//   
//   for (label i=0; i<5; i++)
//   {
//       for (label j=0; j<7; j++)
//       {
//   	Info<< "sourceData_[" << i << "][" << j << "] = " << sourceData_[i][j] << endl;
//       }
//   }

}


Foam::tmp<Foam::volVectorField> Foam::readJBData::assignJBValues
//void Foam::readJBData::assignSourceValues
(
    const fvMesh& mesh
)
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "source",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0),Foam::vector(0,0,0))
        )
    );
    
    volVectorField& source = tsource.ref();
    
    forAll(source,celli)
    {
        //Primeiro, gravo o vetor de coordenadas do ponto central do volume atual
        const Vector<double>& vetorCoordenadas = mesh.C().internalField()[celli];
        //Agora eu procuro no arquivo externo, ja salvo em sourceData_, as coordenadas mais proximas

        double target_x = vetorCoordenadas[0];
        double target_y = vetorCoordenadas[1];
        double target_z = vetorCoordenadas[2];

//       double diferenca_x=1000.0;
//       double diferenca_y=1000.0;
//       double diferenca_z=1000.0;
        double somaDiferenca = 5000.0;

//       label i_x = 0;
//       label i_y = 0;
//       label i_z = 0;
        label a = 0;
        
        for (label i=0; i<nPoints_; i++)
//        forAll(nPoints_,i)
        {
            double temp_x = (target_x - sourceData_[i][1]);
            double temp_y = (target_y - sourceData_[i][2]);
            double temp_z = (target_z - sourceData_[i][3]);

/*            if (temp_x < diferenca_x)
            {
                diferenca_x = temp_x;
            }
            if (temp_y < diferenca_y)
            {
                diferenca_y = temp_y;
            }
            if (temp_z < diferenca_z)
            {
                diferenca_z = temp_z;
            }
*/            
            double tempSoma = pow(pow(temp_x,2) + pow(temp_y,2) + pow(temp_z,2),0.5);
            if (tempSoma < somaDiferenca)
            {
                somaDiferenca = tempSoma;
                a = i;
            }
        }
        source[celli] = Foam::vector(sourceData_[a][4],sourceData_[a][5],sourceData_[a][6]);
    }
//    Info<< "mesh.points() = " << mesh.points()[100] << endl;
//    Info<< "coordenada X " << mesh.points().component(1) << endl;
//    Info<< "coordenada X " << mesh.points()[celula][componente] << endl;

//   Info << "Verificando quais valores esta se colocando para o termo fonte... " << endl;
//   forAll(source,j)
//   //for (label j=0; j<nPoints_; j++)
//   {
//       Info << " -------- Informacoes sobre celula " << j << "--------- " << endl;
//       Info << "Coordenadas:" << endl;
//       Info << "x = " << mesh.C().internalField()[j][0] << endl;
//       Info << "y = " << mesh.C().internalField()[j][1] << endl;
//       Info << "z = " << mesh.C().internalField()[j][2] << endl;
//       Info << "Vetor termo fonte: " << source[j] << endl;
//       Info << " -------------------------------------------" << endl << endl;
//       Info << " -------------------------------------------" << endl << endl;
//       Info << " -------------------------------------------" << endl << endl;
//       
//       Info << "Coordenadas do ponto 0: " << mesh.C().internalField()[0] << endl;
//   }     

    return tsource;
}
