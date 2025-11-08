# Droplet_Mixing

## Description

This repo contains the code, cases files and the data for the simulation of the impact, coalescence and mixing of two (or more) droplets of the same fluid. coalesence can take place in air or on a substrate with the wetting of the droplets able to be modelled using the Kistler contact angle model. 
Impact and coalesence is modelled through a standard VOF approach. Mixing is modelled using a modified diffusion equation that restricts an added scalar $\gamma$ to the droplet phase of the simulation.
Mixing is quantitive assessed through a custome mixing metric which is written out through the simulation in a created postProcessing folder. 

The custom solver and case files for the open-source finite volume library [OpenFOAM 9](https://openfoam.org/), which must first be installed to use this code.

## Problem Parameters
The problem can be characterised by the following physical paramters:

Parameter                             | Symbol        | Case File Location            
--------------------------------------|---------------|-------------------------
Impact velocity                       | $U$           | system/setFieldsDict         
Impacting Droplet radius              | $R$           | system/setFieldsDict
Droplet Viscosity                     | $\mu_{d}$     | constant/transportProperties
Droplet Density                       | $\rho_{d}$    | constant/transportProperties
Gas Viscosity                         | $\mu_{g}$     | constant/transportProperties
Gas Density                           | $\rho_{g}$    | constant/transportProperties
Surface Tension                       | $\sigma$      | constant/transportProperties
Gravity                               | $g$           | constant/g
Receeding Contact Angle               | $\theta_R$    | 0/alpha.drop
Equilibrium Contact Angle             | $\theta_A$    | 0/alpha.drop
Advancing Contact Angle               | $\theta_A$    | 0/alpha.drop
Sesile Droplet Volume                 | $V_s$         | system/setFieldsDict
Sessile Droplet Initial Contact Angle | $\theta_s$    | system/setFieldsDict

Explain how the sessile droplet geometry is set through a psuedo radius etc.

## Installation
This code has been written for OpenFOAM 9, which can be downloaded [here](https://openfoam.org/version/9/)

This repository mirrors the folders in the OpenFOAM installation and therefore can be copied to your user installation.

All customised libraries (solvers and utlities) need compiling with the OpenFOAM command _wmake_ this can also be achieved by cloning the git repo and running the command: _wmake -all_ in the cloned directory
 
## Testing
To test the all custom code has been compiled correctly navigate to a example case file (suggest: X) and run:

```blockMesh``` 

```setFields```

```diffusiveInterFoam```



## Mesh Generation
To create the mesh for the simulation, copy the ```refineSetFields.sh``` script to the example case file and run it using: ```./refineSetFields.sh``` (note: execution privledges may need authorising).
This script will create the inital mesh and set the initial droplet configuration on the generated mesh. 
It is good practice to keep a 0.orig that reflects the state of the ```0``` folder before the ```refineSetFields```
After this the impact simulation can be started with:

```diffusiveInterFoam```


## Static Diffusion
In some cases (such as after the droplets have reached a static equilibrium) it is preferable for simulation speed to stop solving the momentum equations and only solve the diffusion equation with (or without) the residual velocity field. 
This can be done with the ```staticDiffusiveInterFoam``` solver, which only solves the phase restricted diffusion equation over the simulation domain and not the equations of mass and momentum.

In the case the droplets are stationary the ```staticSetFieldsDict``` file can be used in place of the standard ```setFieldsDict``` file to intialise a new simulation with the droplet geometry unalterned but the velocity field in the simulation domain set to zero. 

In this instance then static diffusion of the added scalar takes play only.



## 📁 Repository Structure

solvers/ → source code for custom solvers 

utlities/  → source code for Kistler Model and custom function objects

examples/ → example cases 

data/ → Data from publication

figures/ → Publication plots


## Citation
If you use or build upon this solver, please cite:
```
@software{YatinDarbar_diffusiveInterFoam_2025,
  author = {Yatin Darbar},
  title = {Droplet Mixing -- GitHub Repository},
  year = {2025},
  url = {https://github.com/yatin-darbar/Droplet_Mixing}
}
```
I hope to change this to a publication citation soon!


