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
    Foam::molecularWeightTable

Description
    A table of atomic weights for all the elements

SourceFiles
    molecularWeights.C

\*---------------------------------------------------------------------------*/

#ifndef molecularWeights_H
#define molecularWeights_H

#include "scalar.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class molecularWeights Declaration
\*---------------------------------------------------------------------------*/

class molecularWeightTable
:
    public HashTable<scalar>
{

public:

    // Public types and data

        //- Structure to hold the element name and atomic weight pair
        struct molecularWeight
        {
            char name[5];
            scalar weight;
        };

        static const int nElements = 20;

        //- Static table of the weights of all known elements
        static const molecularWeight molecularWeights[nElements];


    // Constructors

        //- Construct from molecularWeights_
        molecularWeightTable();
};


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

// Atomic weights table for every element in the periodic table
extern molecularWeightTable molecularWeights;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
