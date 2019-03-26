# Multi-species-mass-transport-library-for-OpenFOAM-2.3.X

## Files description

<ul>

<li><p align="justify">"<em>diffusivityModels</em>"</p></li>
    <ul>
        <li><p align="justify">binaryDiffusivityModel -> base class for all binary diffusivity models</p></li>
        <li><p align="justify">constant -> binary diffusivity model with constant diffusion coefficients</p></li>
        <li><p align="justify">Fuller -> binary diffusivity model based on Fuller-Schettler-Giddings correlation</p></li>
        <li><p align="justify">ChapmanEnskog -> binary diffusivity model based of Chapman-Enskog correlation</p></li>
        <li><p align="justify">Wilke -> binary diffusivity model based on Wilke-Lee correlation</p></li>
        <li><p align="justify">Knudsen -> Knudsen diffusivity model</p></li>
        <li><p align="justify">diffusivityModel -> class that collects the binary diffusion coefficients for a set of species</p></li>
        <li><p align="justify">KnudsenDiffusivityModel -> class that collects the Knudsen diffusion coefficients for a set of species</p></li>
    </ul>
    
<li><p align="justify">"<em>multiSpeciesTransportModels</em>"</p></li>
    <ul>
        <li><p align="justify">multiSpeciesTransportModel -> base class for all multiSpecies transport models</p></li>
        <li><p align="justify">Fick -> diffusive mass fluxes are calculated using the common Fick law for multicomponent mixture (D_alpha is function of molar/mass fractions)</p></li>
        <li><p align="justify">FickDilutedMixture -> similar to Fick model, but D_alpha IS NOT function of molar/mass fractions </p></li>
        <li><p align="justify">SchmidtNumber -> diffusive mass fluxes are calculated using the Schmidt number</p></li>
        <li><p align="justify">Lewis -> diffusive mass fluxes are calculated using the Lewis number</p></li>
        <li><p align="justify">Bosanquet -> similar to Fick model, but D_alpha is corrected with Knudsen effects</p></li>
        <li><p align="justify">MaxwellStefan -> diffusive mass fluxes are calculated using the Maxwell-Stefan correlation (j_alpha depend by gradient of all species)</p></li>
    </ul>


<li><p align="justify">"<em>modifiedReactingFoam</em>" in which turbulence is initially generated in the whole domain (usually a box) using Fourier series.</p></li>

</ul>


## Prerequisites

<p align="justify">This library is developed for <strong>OpenFOAM 2.3.X</strong>, and requires the latter version of the Open Source CFD Toolbox. For more information on how to install this version, the reader is referred to <a href="https://sites.google.com/site/foamguides/installation/installing-openfoam-2-3-x">this web page</a>.</p>

## Instructions on program

To compile any solver, go inside the folder of the desired one and execute the command:

    wmake

This will create the application <em>solver</em> in your $FOAM_USER_BIN folder.

All the solvers require adding some additional files:

1. Add thermodynamic pressure level in constant/chemistryProperties
            
        pReff pReff [1 -1 -2 0 0 0 0] 1.01325E+05;

2. A new field-file: pd (dynamic pressure), which is actually the solved pressure variable (p is more or less a dummy field)
