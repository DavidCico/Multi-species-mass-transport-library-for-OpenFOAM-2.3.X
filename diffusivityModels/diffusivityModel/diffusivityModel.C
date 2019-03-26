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

\*---------------------------------------------------------------------------*/

#include "diffusivityModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(diffusivityModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModel::diffusivityModel
(
    const volScalarField& p,
    const volScalarField& T,
    const porosityModelList& pZones,
    const wordList& species
)
:
    dic_
    (
        IOobject
        (
            "transportProperties",
            T.mesh().time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    p_(p),
    
    T_(T),
    
//     eps_(dic_.lookupOrDefault<scalar>("porosity", 1)),
// 
//     tau_(dic_.lookupOrDefault<scalar>("tortuosity", 1)),

    eps_
    (
      IOobject
      (
        "eps",
        T.mesh().time().timeName(),
        T.mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      T.mesh(),
      dimensionedScalar("eps", dimless, 1)
    ),
    
    tau_
    (
      IOobject
      (
        "tau",
        T.mesh().time().timeName(),
        T.mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      T.mesh(),
      dimensionedScalar("tau", dimless, 1)
    ),
    
    pZones_(pZones),
    
    species_(species)
{
    DijModels_.setSize(0.5*species.size()*(species.size()+1));
    Dij_.setSize(DijModels_.size());
    
    forAll(pZones, zoneI) 
    {
        const labelList& cellZoneIds = pZones[zoneI].cellZoneIDs();
        const scalar& porosity = 1;
    forAll(cellZoneIds, zoneJ)
        {
            // const labelList& cells = T.mesh().cellZones()[pZones[zoneI].zoneIds()];
            const labelList& cells = T.mesh().cellZones()[cellZoneIds[zoneJ]];

            forAll(cells, cellI)
            {  
                eps_[cells[cellI]] = porosity;
        }
    }
 }
        
    forAll(species, i)
    {
        for(label j=i; j < species.size(); j++)
        {
            label k = species.size()*i+j-0.5*i*(i+1);
            
            DijModels_.set
            (
                k,
                binaryDiffusivityModel::New
                (
                    species[i],
                    species[j],
                    dic_,
                    p,
                    T
                )
            );

            Dij_.set
            (
                k,
                new volScalarField
                (
                   eps_/tau_*DijModels_[k].D() 
                )
            );   
        }
    } 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diffusivityModel::update()
{    
    for(label i=0; i < species_.size(); i++)
    {
        for(label j=i; j < species_.size(); j++)
        {
            label k = species_.size()*i+j-0.5*i*(i+1);
            Dij_[k] = eps_/tau_*DijModels_[k].D();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
