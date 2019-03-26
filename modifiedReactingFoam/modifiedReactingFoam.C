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

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/
/*Mohsen
#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"*/


#include "fvCFD.H"
#include "fvm.H"
#include "fvOption.H"
#include "fvIOoptionList.H"
#include "meshToMesh.H"
#include "turbulenceModel.H"
#include "psiChemistryCombustion.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"

// ----------------------------- code addition ----------------------------- //
#include "multiSpeciesTransportModel.H"
// ------------------------------------------------------------------------- //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// A.Alexiou 2015, to use OF 2.1 solver style define RHO_EQN_2_1
//#define RHO_EQN_2_1


int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//#   include "readChemistryProperties.H"//////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


    Info<< "Reading chemistry properties\n" << endl;

    IOdictionary chemistryProperties
    (
        IOobject
        (
            "chemistryProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    Switch turbulentReaction(chemistryProperties.lookup("turbulentReaction"));

    dimensionedScalar Cmix("Cmix", dimless, 1.0);

    if (turbulentReaction)
    {
        chemistryProperties.lookup("Cmix") >> Cmix;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


#   include "readGravitationalAcceleration.H"


///////////////////////////////////////////////////////////////////////////////////////////////////////////
#   include "createFields.H"/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


//    //Info<< nl << "Reading thermophysicalProperties" << endl;
//    Info<< nl << "Creating combustion model" << endl;

//    // A.Alexiou 2014
//    // autoPtr<combustionModels:: psiChemistryCombustionModel> combustion // OF 2.1
//    autoPtr<combustionModels::psiCombustionModel> combustion // OF 2.3
//    (
//     // A.Alexiou 2014
//     //combustionModels::psiChemistryCombustionModel::New(mesh) // OF 2.1
//     combustionModels::psiCombustionModel::New(mesh) // OF 2.3
//    );

//    // M. Lindner, A.Alexiou 2014
//    // psiChemistryModel& chemistry = combustion->pChemistry(); // OF 2.1
//    autoPtr<psiChemistryModel> pChemistry // OF 2.3
//    (
//        psiChemistryModel::New(mesh) // OF 2.3
//    );

//    // M. Lindner 2014
//    psiChemistryModel& chemistry = pChemistry(); // OF 2.3

//    // A.Alexiou 2014
//    // hsCombustionThermo& thermo = chemistry.thermo(); // OF 2.1
//    psiReactionThermo& thermo = chemistry.thermo(); // OF 2.3

//    basicMultiComponentMixture& composition = thermo.composition();
//    PtrList<volScalarField>& Y = composition.Y();



//    word inertSpecie(thermo.lookup("inertSpecie"));

//    volScalarField rho
//    (
//        IOobject
//        (
//            "rho",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        thermo.rho()
//    );

//    Info<< "Reading field U\n" << endl;
//    volVectorField U
//    (
//        IOobject
//        (
//            "U",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh
//    );

//    volScalarField& p = thermo.p();
//    const volScalarField& psi = thermo.psi();

//    // A.Alexiou 2014
//    // volScalarField& hs = thermo.hs();
//    volScalarField& hs = thermo.he(); // .hs();

//    // A.Alexiou 2014 - runtime.write() does not work for field T for reason unknown
//    // Instance of an object T didn't help change that
////    volScalarField T
////    (
////        IOobject
////        (
////            "T",
////            runTime.timeName(),
////            mesh,
////            IOobject::MUST_READ,
////            IOobject::AUTO_WRITE
////        ),
////        thermo.T()
////    );

//    const volScalarField& T = thermo.T();

//    #include "compressibleCreatePhi.H"

//    volScalarField kappa
//    (
//        IOobject
//        (
//            "kappa",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("zero", dimless, 0.0)
//    );

//    Info << "Creating turbulence model.\n" << nl;
//    autoPtr<compressible::turbulenceModel> turbulence
//    (
//        compressible::turbulenceModel::New
//        (
//            rho,
//            U,
//            phi,
//            thermo
//        )
//    );

//    // Set the turbulence into the combustion model
//    combustion->setTurbulence(turbulence());

//    /*Info<< "Creating field DpDt\n" << endl;
//    volScalarField DpDt =
//    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);*/

//    Info<< "Creating field dpdt\n" << endl; //Mohsen, may have to use DpDt!
//    volScalarField dpdt("dpdt", fvc::ddt(p));

//    Info<< "Creating field kinetic energy K\n" << endl;
//    volScalarField K("K", 0.5*magSqr(U));

//    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

//    forAll(Y, i)
//    {
//        fields.add(Y[i]);
//    }
//    fields.add(hs);

//    DimensionedField<scalar, volMesh> chemistrySh
//    (
//        IOobject
//        (
//            "chemistry::Sh",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
//    );

//    volScalarField dQ
//    (
//        IOobject
//        (
//            "dQ",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
//    );


//    // ------------------------------------------------------------------------- //
//    autoPtr<multiSpeciesTransportModel> mstm
//    (
//        multiSpeciesTransportModel::New
//        (
//            thermo,
//            turbulence()
//        )
//    );
//    // ------------------------------------------------------------------------- //

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////





#   include "initContinuityErrs.H"
#   include "readTimeControls.H"
#   include "compressibleCourantNo.H"
#   include "setInitialDeltaT.H"

  pimpleControl pimple(mesh);  //Mohsen

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
      //#       include "readPISOControls.H" Mohsen
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"
#include "meshToMesh.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;




///////////////////////////////////////////////////////////////////////////////////////////////////////////
#       include "chemistry.H"////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//        {
//            Info<< "Solving chemistry" << endl;

////                    A.Alexiou 2015
////                    chemistry.solve
////                    (
////                        runTime.value() - runTime.deltaT().value(),
////                        runTime.deltaT().value()
////                    );

//            chemistry.solve
//            (
//                runTime.deltaT().value()
//            );


//            // turbulent time scale
//            if (turbulentReaction)
//            {
//                volScalarField tk =
//                        Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
//                volScalarField tc = chemistry.tc();

//                // Chalmers PaSR model
//                kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
//            }
//            else
//            {
//                kappa = 1.0;
//            }

//            chemistrySh = kappa*chemistry.Sh()();
//        }


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

        // M.Lindner 2014, A.Alexiou 2015
        fv::IOoptionList fvOptions(mesh); // for compatibility with OF 2.3, since rhoEqn.H changed


#ifndef RHO_EQN_2_1

        // A.Alexiou 2015
#include "rhoEqn.H" // OF header/source

#else

// A.Alexiou 2015
// original rhoEqn.H from OF 2.1, works fine with OF 2.3, but then rhoEqn.H has to be
// commented and 'fv::IOoptionList fvOptions(mesh);' is not used

        {
            solve(fvm::ddt(rho) + fvc::div(phi));
        }

#endif

// A.Alexiou - for debugging
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//#include "rhoEqn.H"//////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//      A.Alexiou 2015, for comparison to OF 2.1 rhoEqn.H included and commented

//        {
//            fvScalarMatrix rhoEqn
//            (
//                fvm::ddt(rho)
//              + fvc::div(phi)
//              ==
//                fvOptions(rho)
//            );

//            fvOptions.constrain(rhoEqn);

//            rhoEqn.solve();

//            fvOptions.correct(rho);
//        }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


        //for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)  Mohsen
        while (pimple.loop()) //Mohsen
        {


///////////////////////////////////////////////////////////////////////////////////////////////////////////
#           include "UEqn.H"/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


//            fvVectorMatrix UEqn
//                    (
//                        fvm::ddt(rho, U)
//                        + fvm::div(phi, U)
//                        + turbulence->divDevRhoReff(U)
//                        ==
//                        rho*g
//                        );

//            UEqn.relax();

//            /* if (momentumPredictor)  Mohsen
//            {
//                solve(UEqn == -fvc::grad(p));
//            }*/
//            if (pimple.momentumPredictor()) //Mohsen
//            {
//                solve(UEqn == -fvc::grad(p));
//                K = 0.5*magSqr(U);
//            }


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////
#           include "YEqn.H"/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//            tmp<fv::convectionScheme<scalar> > mvConvection
//                    (
//                        fv::convectionScheme<scalar>::New
//                        (
//                            mesh,
//                            fields,
//                            phi,
//                            mesh.divScheme("div(phi,Yi_h)")
//                            )
//                        );


//            // A.Alexiou 2015
//            // mstm().correct(kappa, chemistry, fields); // OF 2.1
//            // ------------------------------------------------------------------------- //
//            // A.Alexiou 2015
//            mstm().correct(Y, kappa, chemistry, fields); // OF 2.3
//            // ------------------------------------------------------------------------- //

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////
#           include "hsEqn.H"////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



//            {
//                fvScalarMatrix hsEqn
//                        (
//                            fvm::ddt(rho, hs)
//                            + mvConvection->fvmDiv(phi, hs)
//                            - fvm::laplacian(turbulence->alphaEff(), hs)
//                            + mstm().multiSpeciesHeatSource()
//                            ==
//                            dpdt
//                            // ------------------------------------------------------------------------- //
//                            + chemistrySh
//                            // ------------------------------------------------------------------------- //
//                            );

//                hsEqn.relax();
//                hsEqn.solve();

//                thermo.correct();

//                // A.Alexiou 2015 - Debugging
//                //Info << "Check T" << thermo.T();

//                Info<< "T gas min/max   = " << min(T).value() << ", "
//                    << max(T).value() << endl;
//            }



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



            // --- PISO loop
            //for (int corr=1; corr<=nCorr; corr++) Mohsen
            while (pimple.correct())
            {



// A.Alexiou - Changed pEqn from OF 2.1 pEqn.H header file to OF 2.3 pEqn.H
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#               include "pEqn.H"/////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//                rho = thermo.rho();

//                volScalarField rAU(1.0/UEqn.A());
//                surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));  // OF 2.3

//                volVectorField HbyA("HbyA",U);  // OF 2.3
//                HbyA = rAU*UEqn.H(); // OF 2.3

//                if (pimple.transonic())
//                {

//                    //               // 2.1
//                    //               surfaceScalarField phid
//                    //               (
//                    //                   "phid",
//                    //                   fvc::interpolate(psi)
//                    //                  *(
//                    //                       (fvc::interpolate(U) & mesh.Sf())
//                    //                     + fvc::ddtPhiCorr(rAU, rho, U, phi)
//                    //                   )
//                    //               );


//                    // 2.3
//                    surfaceScalarField phid
//                            (
//                                "phid",
//                                fvc::interpolate(psi)
//                                *(
//                                    (fvc::interpolate(rho*HbyA) & mesh.Sf())
//                                    + rhorAUf*fvc::ddtCorr(rho, U, phi)
//                                    )/fvc::interpolate(rho)
//                                );

//                    fvOptions.makeRelative(fvc::interpolate(psi), phid);


//                    //               // 2.1
//                    //               while (pimple.correctNonOrthogonal())
//                    //               {
//                    //                   fvScalarMatrix pEqn
//                    //                   (
//                    //                       fvm::ddt(psi, p)
//                    //                     + fvm::div(phid, p)
//                    //                     - fvm::laplacian(rho*rAU, p)
//                    //                   );

//                    //               pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

//                    //               if (pimple.finalNonOrthogonalIter())
//                    //                   {
//                    //                       phi == pEqn.flux();
//                    //                   }
//                    //               }


//                    while (pimple.correctNonOrthogonal())
//                    {
//                        fvScalarMatrix pEqn
//                                (
//                                    fvm::ddt(psi, p)
//                                    + fvm::div(phid, p)
//                                    - fvm::laplacian(rho*rAU, p)
//                                    ==
//                                    fvOptions(psi, p, rho.name())
//                                    );

//                        fvOptions.constrain(pEqn);

//                        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

//                        if (pimple.finalNonOrthogonalIter())
//                        {
//                            phi == pEqn.flux();
//                        }
//                    }


//                }
//                else
//                {

//                    //               // 2.1
//                    //               phi =
//                    //                   fvc::interpolate(rho)
//                    //                  *(
//                    //                       (fvc::interpolate(U) & mesh.Sf())
//                    //                     + fvc::ddtPhiCorr(rAU, rho, U, phi)
//                    //                   );


//                    // 2.3
//                    surfaceScalarField phiHbyA
//                            (
//                                "phiHbyA",
//                                (
//                                    (fvc::interpolate(rho*HbyA) & mesh.Sf())
//                                    + rhorAUf*fvc::ddtCorr(rho, U, phi)
//                                    )
//                                );

//                    fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);


//                    //               // 2.1
//                    //               while (pimple.correctNonOrthogonal())
//                    //               {
//                    //                   fvScalarMatrix pEqn
//                    //                   (
//                    //                       fvm::ddt(psi, p)
//                    //                     + fvc::div(phi)
//                    //                     - fvm::laplacian(rho*rAU, p)
//                    //                   );

//                    //               pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

//                    //               if (pimple.finalNonOrthogonalIter())
//                    //                   {
//                    //                       phi += pEqn.flux();
//                    //                   }
//                    //               }

//                    // 2.3
//                    while (pimple.correctNonOrthogonal())
//                    {
//                        fvScalarMatrix pEqn
//                                (
//                                    fvm::ddt(psi, p)
//                                    + fvc::div(phiHbyA)
//                                    - fvm::laplacian(rho*rAU, p)
//                                    ==
//                                    fvOptions(psi, p, rho.name())
//                                    );

//                        fvOptions.constrain(pEqn);

//                        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

//                        if (pimple.finalNonOrthogonalIter())
//                        {
//                            phi = phiHbyA + pEqn.flux();
//                        }
//                    }



//                }


//#ifndef RHO_EQN_2_1

//#include "rhoEqn.H"

//#else

//                // A.Alexiou 2015
//                // original rhoEqn.H from OF 2.1, works fine with OF 2.3, but then rhoEqn.H has to be
//                // commented and 'fv::IOoptionList fvOptions(mesh);' is not used
//                {
//                    solve(fvm::ddt(rho) + fvc::div(phi));
//                }

//#endif

//                //           // 2.1
//                //           #include "compressibleContinuityErrs.H"

//                //           U -= rAU*fvc::grad(p);
//                //           U.correctBoundaryConditions();

//                //           K = 0.5*magSqr(U);

//                //           //DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
//                //           dpdt = fvc::ddt(p);


//#include "compressibleContinuityErrs.H"

//                U = HbyA - rAU*fvc::grad(p);
//                U.correctBoundaryConditions();
//                fvOptions.correct(U);
//                K = 0.5*magSqr(U);

//                if (thermo.dpdt())
//                {
//                    dpdt = fvc::ddt(p);
//                }


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////




            }
            /*} Mohsen

         turbulence->correct();*/

            // A.Alexiou 2015
//            if (runTime.write())
//            {
//                chemistry.dQ()().write();
//            }
            if (pimple.turbCorr())  //Mohsen
            {
                turbulence->correct();
            }
        }

        // A.Alexiou 2015
        //runTime.write();
        if (runTime.write())
        {
            // A.Alexiou 2015
            thermo.T().write();
            chemistry.dQ()().write();
            // A.Alexiou 2015
            forAll(composition.Y(), i)
            {
                composition.Y()[i].write();
            }
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
