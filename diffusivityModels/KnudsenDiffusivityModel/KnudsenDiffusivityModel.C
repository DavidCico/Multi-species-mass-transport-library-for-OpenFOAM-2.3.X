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

#include "KnudsenDiffusivityModel.H"

namespace Foam
{
    defineTypeNameAndDebug(KnudsenDiffusivityModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KnudsenDiffusivityModel::KnudsenDiffusivityModel
(
    const volScalarField& T,
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
    
    T_(T),
    
    eps_(dic_.lookupOrDefault<scalar>("porosity", 1)),

    tau_(dic_.lookupOrDefault<scalar>("tortuosity", 1)),
    
    species_(species)
{
    DKModels_.setSize(species.size());
    DK_.setSize(DKModels_.size());
    
    forAll(species, i)
    {
        DKModels_.set
        (
            i,
            new Knudsen
            (
                species[i],
                dic_,
                T
            )
        );

        DK_.set
        (
            i,
            new volScalarField
            (
                eps_/tau_*DKModels_[i].DK() 
            )
        );
    } 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::KnudsenDiffusivityModel::update()
{    
    forAll(species_, i)
    {
        DK_[i] = eps_/tau_*DKModels_[i].DK();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
