/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::multiSpeciesTransportModel

Description
    Abstract base class for multispecies mass transport model models.

SourceFiles
    multiSpeciesTransportModelI.H
    multiSpeciesTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiSpeciesTransportModel_H
#define multiSpeciesTransportModel_H

#include "fvc.H"

#include "psiChemistryModel.H"
#include "turbulenceModel.H"
#include "diffusivityModel.H"
#include "multiComponentMixture.H"
#include "porosityModelList.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                    Class multiSpeciesTransportModel Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesTransportModel
:

    public IOdictionary
{
    
protected:  
  
    // Protected data

        //- Mesh
        const fvMesh& mesh_;

        // Porous zones
        // M.Lindner, A.Alexiou 16.12.2014
        // porousZones pZones_;
        porosityModelList pZones_ ;

        // Reference to the thermodynamic package
        // M.Lindner, A.Alexiou 16.12.2014
        // hsCombustionThermo& thermo_;
        psiReactionThermo& thermo_;

        //- Index of onert species
        const scalar inertIndex_;

        //- Reference to the turbulence model
        const compressible::turbulenceModel& turbulence_;

        //- Molar fractions [-]
        PtrList<volScalarField> x_;

        //- Species mass fluxes [kg/s]
        PtrList<surfaceScalarField> n_;

        //- Species rates [kg/m^3/s]
        PtrList<volScalarField> Sy_;

        //- Diffusivity models 
        autoPtr<diffusivityModel> DijModel_;

        //- Multispecies diffusion coefficients
        PtrList<volScalarField> D_;

        //- Costruct molar fraction fields
        void updateMolarFractions();


private:
  
    // Private member functions
      
        //- Disallow copy construct
        multiSpeciesTransportModel(const multiSpeciesTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesTransportModel&);

        tmp<Foam::volScalarField> RR
        (
            const word Yname,
            const psiChemistryModel& chemistry,
            const label i
        ) const;

//        Foam::multiSpeciesTransportModel::RR
//        (
//            const word Yname,
//            const psiChemistryModel& chemistry,
//            const label i
//        ) const


public:
 
    //- Runtime type information
    TypeName("multiSpeciesTransportModel");

      
    // Declare run-time constructor selection table
       
#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        multiSpeciesTransportModel,
        fvMesh,
        (
            psiReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        ),
        (thermo, turbulence)
    );
#endif
  

    // Constructors	  

        //- Construct from components
        multiSpeciesTransportModel
        (
            psiReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    // Selectors	  
  
        //- Return a reference to the selected multispecies transport model

        static autoPtr<multiSpeciesTransportModel> New
        (
            psiReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    //- Destructor

        virtual ~multiSpeciesTransportModel()
        {}


    // Member Functions

        // Access   
           
            //- Return the const binary diffusion coefficient field
            inline const volScalarField& Dij(label i, label j) const;

            //- Return the porous zones
            inline porosityModelList& pZones();
            
            //- Return the table of species
            inline const speciesTable& species() const;

            //- Return the const total mass flux field for a species
            inline const surfaceScalarField& n(const label i) const;

            //- Return the diffusive mass flux field for a specie given by index
            surfaceScalarField j(const label i);

            //- Molar weight for a specie given by index
            virtual scalar W(label i) const = 0;

            //- Specific enthalpy for a specie given by index and temperature
            virtual scalar hs(label i, label T) const = 0;


        //- Sorce term in hentalpy equation due to diffusive fluxes
        tmp<volScalarField> multiSpeciesHeatSource();

        //- Solve the mass transport equantions without reactions
        virtual scalar correct
        (
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        ) = 0;

        //- Solve the mass transport equantions with reactions
        scalar correct
        (
            const PtrList<volScalarField>& Y,
            const volScalarField& kappa,
            const psiChemistryModel& chemistry,
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        );

        //- Read transportProperties dictionary
        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiSpeciesTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
