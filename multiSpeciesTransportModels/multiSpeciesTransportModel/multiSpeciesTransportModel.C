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

#include "multiSpeciesTransportModel.H"
#include "dimensionedConstants.H"
#include "constants.H" //Mohsen
#include "psiChemistryCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(multiSpeciesTransportModel, 0);
    defineRunTimeSelectionTable(multiSpeciesTransportModel, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiSpeciesTransportModel::updateMolarFractions()
{
  //const scalar R = dimensionedConstant("R", 8314.51); //Mohsen
  const scalar R = constant::physicoChemical::R.value()*1000;

    forAll(thermo_.composition().Y(), i)
    {
        const scalarField& psiCells = thermo_.psi().internalField();
        const scalarField& TCells = thermo_.T().internalField();
        const scalarField& yCells =
            thermo_.composition().Y(i).internalField();

        scalarField& xCells = x_[i].internalField();

        forAll(xCells, cellI)
        {
           xCells[cellI] = yCells[cellI] * psiCells[cellI] * TCells[cellI]
                * R / W(i);
        }
        forAll(x_[i].boundaryField(), patchI)
        {
            const fvPatchScalarField& ppsi =
                thermo_.psi().boundaryField()[patchI];

            const fvPatchScalarField& pT =
                thermo_.T().boundaryField()[patchI];

            const fvPatchScalarField& py =
                thermo_.composition().Y(i).boundaryField()[patchI];

            fvPatchScalarField& px = x_[i].boundaryField()[patchI];

            forAll(px, faceI)
            {
                px[faceI] = py[faceI] * ppsi[faceI] * pT[faceI]
                    * R / W(i);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSpeciesTransportModel::multiSpeciesTransportModel
(
    psiReactionThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE

        )
    ),

    mesh_(thermo.T().mesh()),

    // M.Lindner, A.Alexiou 16.12.2014
    pZones_
    (
        thermo.T().mesh(),
        IOdictionary
        (
            IOobject
            (
                "porousZones",
                thermo.T().mesh().time().constant(),
                thermo.T().mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),

    thermo_(thermo),

    // The last specie in the list is the default inertspecie
    inertIndex_(species().size()-1),

    turbulence_(turbulence)
{
    // Construct the diffusivity model
    DijModel_.set(new diffusivityModel(thermo.p(), thermo.T(), pZones_, species()));

    x_.setSize(species().size());
    n_.setSize(species().size());
    Sy_.setSize(species().size());

    forAll(thermo.composition().Y(), i)
    {
        x_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "x_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimless
            )
        );

        n_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "n_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("n", dimMass/dimTime, 0.0)
            )
        );

        Sy_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "Sy_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("Sy", dimMass/dimTime/dimVolume, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::surfaceScalarField Foam::multiSpeciesTransportModel::j
(
    const Foam::label i
)
{
    return n(i) -
        fvc::interpolate(thermo_.composition().Y(i)) * turbulence_.phi();
}


// A.Alexiou 2015
// Former function correct running with OF 2.1
//Foam::scalar Foam::multiSpeciesTransportModel::correct
//(
//    volScalarField& kappa,
//    const psiChemistryModel& chemistry,
//    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
//)
//{

//    forAll(Sy_, i)
//    {
//          Sy_[i] = kappa*chemistry.RR(i);
            //Sy_[i] = chemistry.RR(i);//Mohsen
//    }

//    return correct(fields);
//}


// A.Alexiou 2015
// Modified function correct for OF 2.3
Foam::scalar Foam::multiSpeciesTransportModel::correct
(
    const PtrList<volScalarField>& Y,
    const volScalarField& kappa,
    const psiChemistryModel& chemistry,
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{

    forAll(Sy_, i)
    {
          // A.Alexiou 2015
          // Sy_[i] = kappa*chemistry.RR(i); // OF 2.1
          Sy_[i] = kappa*RR(Y[i].name(), chemistry,i); // OF 2.3
    }

    return correct(fields);
}



// A.Alexiou 2015
// Added from ODEChemistryModel OF 2.1 and modified for OF 2.3
//template<class CompType, class ThermoType>
inline Foam::tmp<Foam::volScalarField>
//Foam::ODEChemistryModel<CompType, ThermoType>::RR
Foam::multiSpeciesTransportModel::RR
(
    const word Yname,
    const psiChemistryModel& chemistry,
    const label i
) const
{

    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR(" + Yname + ')',
                chemistry.time().timeName(),
                chemistry.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            chemistry.mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

//    // M. Yusufi, A.Alexiou 2015
//    volScalarField& ttRR = tRR();

    if (chemistry.chemistry())
    {
        // A.Alexiou 2015
        tRR().internalField() = chemistry.RR(i);
        tRR().correctBoundaryConditions();
//        ttRR.internalField() = chemistry.RR(i);
//        ttRR.correctBoundaryConditions();
    }
    return tRR;
}



Foam::tmp<Foam::volScalarField>
Foam::multiSpeciesTransportModel::multiSpeciesHeatSource()
{
    tmp<volScalarField> tMsht
    (
        new volScalarField
        (
            IOobject
            (
                "multiSpeciesHeatSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("mths", dimensionSet(1, -1, -3, 0, 0), 0)
        )
    );

    volScalarField& msht = tMsht();

    volScalarField hsN("hsN",thermo_.he());

    forAll(hsN, celli)
    {
        hsN[celli] = hs(inertIndex_,thermo_.T()[celli]);
    }

    forAll(hsN.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT =
            thermo_.T().boundaryField()[patchi];

        fvPatchScalarField& pHsN = hsN.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pHsN[facei] = hs(inertIndex_,pT[facei]);
        }
    }

    forAll(species(), i)
    {
        if (i != inertIndex_)
        {

        volScalarField hsi("hsi",thermo_.he());

        forAll(hsi, celli)
        {
            hsi[celli] = hs(i,thermo_.T()[celli]);
        }

        forAll(hsi.boundaryField(), patchi)
        {
            const fvPatchScalarField& pT =
                thermo_.T().boundaryField()[patchi];

            fvPatchScalarField& pHsi = hsi.boundaryField()[patchi];

            forAll(pT, facei)
            {
                pHsi[facei] = hs(i,pT[facei]);
            }
        }

//         volScalarField alphaMinusRhoDHsi = (turbulence_.alphaEff() - D_[i]) * (hsi - hsN);
        volScalarField alphaMinusRhoDHsi = (turbulence_.alphaEff()) * (hsi - hsN);

        volScalarField laplacianAlpHsi =
            fvc::laplacian(alphaMinusRhoDHsi, thermo_.composition().Y(i));

        volScalarField deltahs = hsi-hsN;
        laplacianAlpHsi += fvc::div(j(i),deltahs, "div(ji,hi)");

        forAll(msht, celli)
        {
            msht[celli] += laplacianAlpHsi[celli];
        }

        forAll(msht.boundaryField(), patchi)
        {
            const fvPatchScalarField& plaplacianAlpHsi =
                laplacianAlpHsi.boundaryField()[patchi];

            fvPatchScalarField& pmsht = msht.boundaryField()[patchi];

            forAll(pmsht, facei)
            {
                pmsht[facei] += plaplacianAlpHsi[facei];
            }
        }

        }
    }

    return tMsht;
}


bool Foam::multiSpeciesTransportModel::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
