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

#include "Knudsen.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  
    defineTypeNameAndDebug(Knudsen, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Knudsen::Knudsen
(
    const word& name,
    const dictionary& dic,
    const volScalarField& T
)
:
    dic_(dic),
    
    name_(name),

    T_(T)
{  
    d = dimensionedScalar(dic_.lookup("poreDiametre")).value();
  
    //W = molecularWeights[name];
    const scalar& W = readScalar(dic.subDict("molarWeight").lookup(name));
    
    RR = 8314.51;
    PI = 3.14159;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::Knudsen::DK
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tDk(new scalarField(T.size()));
    scalarField& dk = tDk();

    forAll(T, facei)
    {
        dk[facei] = d / 3 * sqrt((8*RR*T[facei])/(PI*W));
    }

    return tDk;
}


Foam::tmp<Foam::volScalarField> Foam::Knudsen::DK() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tDk
    (
        new volScalarField
        (
            IOobject
            (
                "Dk_" + name_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -1, 0, 0)
        )
    );

    volScalarField& dk = tDk();

    forAll(this->T_, celli)
    {
	dk[celli] = d / 3 * sqrt((8*RR*T_[celli])/(PI*W));
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pDk = dk.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pDk[facei] = d / 3 * sqrt((8*RR*T_[facei])/(PI*W));
        }
    }

    return tDk;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
