# Multi-species-mass-transport-library-for-OpenFOAM-2.3.X

<p align="justify">This repository contains different solvers for OpenFOAM 2.3.X. All of the different solvers, are used for combustion and use a <a href="https://www.sharcnet.ca/Software/Ansys/17.0/en-us/help/cfx_thry/i1309364.html"> finite rate chemistry approach </a>, in which the different species and elementary chemical reactions are modelled following the Arrhenius law.</p>

<p align="justify">The different solvers are based on a <a href="https://ccse.lbl.gov/Research/LowMach/lowMach.html"> low Mach number</a> assumption, where the pressure is decoupled from the thermodynamic properties of the fluid.</p>

Duwig, Christophe, et al. Large Eddy Simulations of a piloted lean premix jet flame using finite-rate chemistry. <em>Combustion Theory and Modelling</em>, 15(4):537-568, 2011.

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
<li><p align="justify">"<em>dynreactingLMFoam</em>" which accomodates of a moving mesh.</p></li>
<li><p align="justify">"<em>reactingBoxturbFoam</em>" in which turbulence is initially generated in the whole domain (usually a box) using Fourier series.</p></li>

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
