# Droplet_Mixing

## Description

This repo contains the code, cases files and the data for the simulation of the impact, coalescence and mixing of two (or more) droplets of the same fluid. coalesence can take place in air or on a substrate with the wetting of the droplets able to be modelled using the Kistler contact angle model. 

The custom solver and case files for the open-source finite volume library [OpenFOAM 9](https://openfoam.org/), which must first be installed to use this code.

## Problem Parameters
The problem can be characterised by the following physical paramters:

Parameter                             | Symbol        | Case File Location
---------------------------------------------------------------------------------------
Impact velocity                       | $U$           | system/setFieldsDict
Impacting Droplet radius              | $R$           | system/setFieldsDict
Droplet Viscosity                     | $\mu_{d}      | constant/transportProperties
Droplet Density                       | $\rho_{d}     | constant/transportProperties
Gas Viscosity                         | $\mu_{g}      | constant/transportProperties
Gas Density                           | $\rho_{g}     | constant/transportProperties
Surface Tension                       | $\sigma$      | constant/transportProperties
Gravity                               | $g$           | constant/g
Receeding Contact Angle               | $\theta_R$    | 0/alpha.drop
Equilibrium Contact Angle             | $\theta_A$    | 0/alpha.drop
Advancing Contact Angle               | $\theta_A$    | 0/alpha.drop
Sesile Droplet Volume                 | $V_s$         | system/setFieldsDict
Sessile Droplet Initial Contact Angle | $\theta_s$    | system/setFieldsDict

Explain how the sessile droplet geometry is set through a psuedo radius etc.

