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

#include "Fuller.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusivityModels
    {
        defineTypeNameAndDebug(Fuller, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusivityModel,
            Fuller, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModels::Fuller::Fuller
(
    const word& name1,
    const word& name2,
    const dictionary& dic,
    const volScalarField& p,
    const volScalarField& T
)
:
    binaryDiffusivityModel(name1, name2, dic, p, T)
{   
  //const scalar& W1 = molecularWeights[name1];
  //const scalar& W2 = molecularWeights[name2];
  const scalar& W1 = readScalar(dic.subDict("molarWeight").lookup(name1));
  const scalar& W2 = readScalar(dic.subDict("molarWeight").lookup(name2));

    const scalar& V1 = readScalar(dic.subDict("Fuller").lookup(name1));
    const scalar& V2 = readScalar(dic.subDict("Fuller").lookup(name2));

    W12 = (W1 * W2) / (W1 + W2);
    V12 = sqr( pow( V1, (1.0/3.0) ) + pow( V2, (1.0/3.0)) );
    phi =  101325;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModels::Fuller::D
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& d = tD();

    forAll(T, facei)
    {
        d[facei] =
            1.011e-7 * phi * pow(T[facei],(1.75)) / sqrt(W12) / p[facei] / V12;
    }

    return tD;
}


Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusivityModels::Fuller::D() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                "D_" + name1_ + "_" + name2_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -1, 0, 0)
        )
    );

    volScalarField& d = tD();

    forAll(this->T_, celli)
    {
        d[celli] =
            1.011e-7 * phi * pow(T_[celli],(1.75))
            / sqrt(W12) / p_[celli] / V12;
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pD[facei] =
                1.011e-7 * phi * pow(pT[facei],(1.75))
                / sqrt(W12) / pp[facei] / V12;
        }
    }

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
