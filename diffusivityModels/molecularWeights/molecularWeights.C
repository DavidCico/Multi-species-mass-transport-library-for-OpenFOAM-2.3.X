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

#include "molecularWeights.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::molecularWeightTable::molecularWeight
Foam::molecularWeightTable::molecularWeights[molecularWeightTable::nElements] = 
{
    {"H",     1.00797},
    {"O",    15.99940},
    {"H2",    2.01594},
    {"O2",   31.99880},    
    {"H2O",  18.01534},
    {"H2O2", 34.01474},
    {"OH",   17.00737},    
    {"HO2",  33.00677},   
    {"CO2",  44.00995},    
    {"CH4",  16.04303},  
    {"N2",   28.01340},   
    {"AR",   39.948},
    {"",      0.0},
    {"",      0.0},
    {"",      0.0}
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molecularWeightTable::molecularWeightTable()
{
    for (int i=0; i<nElements; i++)
    {
        insert(word(molecularWeights[i].name), molecularWeights[i].weight);
    }
}


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

Foam::molecularWeightTable Foam::molecularWeights;


// ************************************************************************* //
