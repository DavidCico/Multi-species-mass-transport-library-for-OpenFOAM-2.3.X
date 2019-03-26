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
    Foam::MaxwellStefan

Description
    Maxwell-Stefan multispecies mass transport model.
  
SourceFiles
    MaxwellStefan.C

\*---------------------------------------------------------------------------*/

#ifndef MaxwellStefan_H
#define MaxwellStefan_H

#include "multiSpeciesTransportModel.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                         Class MaxwellStefan Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class MaxwellStefan
:
    public multiSpeciesTransportModel
{     
  
protected:
  
    // Protected data

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;

        //- Inverse of multispecies diffusion coefficients (for boundary condition)
        PtrList<volScalarField> G_;


    // Protected member functions
    
        //- Update the coefficients of the species equations
        void updateCoefficients(); 
    

public:    

    //- Runtime type information
    TypeName("MaxwellStefan");

    
    // Constructors

        //- Construct from components
        MaxwellStefan
        (
            // A.Alexiou 2014
            // hsCombustionThermo& thermo,
            psiReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    //- Destructor
    ~MaxwellStefan()
    {}


    // Member Functions

        //- Molar weight for a specie given by index
        scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }
    
        //- Specific enthalpy for a specie given by index and temperature
        scalar hs(label i, label T) const
        {
            // A.Alexiou 2014
            // return speciesThermo_[i].Hs(T)
            return speciesThermo_[i].hs(i,T);
        }
        
        inline volScalarField& D(label i, label j)
        {
            label k = (species().size()-1)*i+j;
            return D_[k];
        }

        inline volScalarField& G(label i, label j)
        {
            label k = (species().size()-1)*i+j;
            return G_[k];
        }
        
        //- Solve the mass transport equantions without reactions
        virtual scalar correct
        (
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        );
                
        //- Read transportProperties dictionary
        virtual bool read();
  };
  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "MaxwellStefan.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
