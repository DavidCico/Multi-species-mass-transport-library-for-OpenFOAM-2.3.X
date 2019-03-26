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
    Foam::Knudsen

Description
    Knudsen diffusivity model.
    
    Knudsen diffusivity coefficient is obtained by
    @f[
        D_{ik} = r \cdot \frac{2}{3} \cdot \sqrt{\frac{8 R T}{\pi M_i}}
    @f]

SourceFiles
    Knudsen.C

\*---------------------------------------------------------------------------*/

#ifndef Knudsen_H
#define Knudsen_H

#include "specie.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "volFields.H"
#include "molecularWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Knudsen Declaration
\*---------------------------------------------------------------------------*/

class Knudsen
{
    // Private data
    
        //- Transport propertied dictionary
        const dictionary& dic_;

        //- Specie name
        const word& name_;

        // Temperature [K]
        const volScalarField& T_;

        //- Perfect gas constant value (J/kmol/K)
        scalar RR, PI;

        //- Specie molar weight
        scalar W;

        //- Mean pore diametre
        scalar d;
    

public:

    //- Runtime type information
    TypeName("Knudsen");


    // Constructors

        //- construct from components
        Knudsen
        (
            const word& name,
            const dictionary& transportDic,
            const volScalarField& T
        );


    // Destructor

        virtual ~Knudsen()
        {}


    // Member Functions
        
        //- Knudsen diffusion coefficient for patch [m^2/s]
        virtual tmp<scalarField> DK
        (
            const scalarField& T,
            const label patchi
        ) const;

        //- Knudsen diffusion coefficient [m^2/s]
        virtual tmp<volScalarField> DK() const;                
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
