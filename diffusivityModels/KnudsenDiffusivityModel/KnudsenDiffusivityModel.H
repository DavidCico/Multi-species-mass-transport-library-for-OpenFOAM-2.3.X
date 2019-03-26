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
    Foam::KnudsenDiffusivityModel

Description
    Diffusivity model for multispecies Knudsen diffusion coefficients 

SourceFiles
    KnudsenDiffusivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef KnudsenDiffusivityModel_H
#define KnudsenDiffusivityModel_H

#include "Knudsen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class KnudsenDiffusivityModel Declaration
\*---------------------------------------------------------------------------*/

class KnudsenDiffusivityModel
{

protected:

    // Protected data
    
        //- Transport propertied dictionary
        const IOdictionary dic_;

        // Temperature [K]
        const volScalarField& T_;

        //- Porosity (if present)
        scalar eps_;

        //- Tortuosity (if present)
        scalar tau_;

        wordList species_;

        //- Knudsen diffusivity models 
        PtrList<Knudsen> DKModels_;   

        //- Knudsen diffusion coefficients
        PtrList<volScalarField> DK_;  


    // Private Member Functions

        //- Disallow copy construct
        KnudsenDiffusivityModel(const KnudsenDiffusivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const KnudsenDiffusivityModel&);


public:

    //- Runtime type information
    TypeName("KnudsenDiffusivityModel");

    // Constructors

        //- Construct from components
        KnudsenDiffusivityModel
        (
            const volScalarField& T,
            const wordList& species
        );


    // Destructor

        virtual ~KnudsenDiffusivityModel()
        {}


        // Member Functions

        inline const volScalarField& D(label i) const
        {
            return DK_[i];
        }

    
        //- Update the diffusion coefficients
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
