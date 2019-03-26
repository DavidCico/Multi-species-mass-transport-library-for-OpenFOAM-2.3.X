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

#include "Wilke.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusivityModels
    {
        defineTypeNameAndDebug(Wilke, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusivityModel,
            Wilke, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModels::Wilke::Wilke
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
    
    A = 1.06036;    B = 0.15610;    C = 0.19300;    DD = 0.47635;
    E = 1.03587;    F = 1.52996;    G = 1.76474;    H = 3.89411;
    
    //const scalar& W1 = molecularWeights[name1];
    //const scalar& W2 = molecularWeights[name2];

    const scalar& W1 = readScalar(dic.subDict("molarWeight").lookup(name1));
    const scalar& W2 = readScalar(dic.subDict("molarWeight").lookup(name2));

    const scalar& epsLJ1 = readScalar(dic.subDict("epsilonLJ").lookup(name1));
    const scalar& epsLJ2 = readScalar(dic.subDict("epsilonLJ").lookup(name2));
    
    const scalar& sigma1 =
        readScalar(dic.subDict("collisionalDiametre").lookup(name1));
    const scalar& sigma2 =
        readScalar(dic.subDict("collisionalDiametre").lookup(name2));  

    sigma_ij = (sigma1 + sigma2) / 2;
    W12 = (W1 * W2) / (W1 + W2);
    phi =  101325; 
    E_ij = sqrt(epsLJ1*epsLJ2);

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModels::Wilke::D
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
        scalar T_N = T[facei] / E_ij;
        scalar Omega =
            A/pow(T_N,B) + C/exp(DD*T_N) + E/exp(F*T_N) + G/exp(H*T_N);
        d[facei] =
            1.858e-7 * phi * sqrt(pow(T[facei],3) / W12)
            / p[facei] / sqr(sigma_ij) / Omega;
    }

    return tD;
}


Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusivityModels::Wilke::D() const
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
        scalar T_N = T_[celli] / E_ij;
        scalar Omega =
            A/pow(T_N,B) + C/exp(DD*T_N) + E/exp(F*T_N) + G/exp(H*T_N);
        d[celli] =
            1.858e-7 * phi * sqrt(pow(T_[celli],3) / W12)
            / p_[celli] / sqr(sigma_ij) / Omega;
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryField()[patchi];

        forAll(pT, facei)
        {
            scalar T_N = pT[facei] / E_ij;
            scalar Omega =
                A/pow(T_N,B) + C/exp(DD*T_N) + E/exp(F*T_N) + G/exp(H*T_N);
            pD[facei] =
                1.858e-7 * phi * sqrt(pow(pT[facei],3) / W12)
                / pp[facei] / sqr(sigma_ij) / Omega;
        }
    }

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
