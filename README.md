# 💧 Droplet_Mixing

## 🧠 Overview

This repository contains the **custom OpenFOAM solvers, case files, and data** used to simulate the **impact, coalescence, and mixing** of two (_or more_) droplets composed of the same fluid.

Droplet interactions can occur either **in air** or **on a solid substrate**, with wetting behavior modeled using the **Kistler contact angle model**.

- **Impact** and **coalescence** are modeled using a **Volume of Fluid** approach.
- **Mixing** is captured via a **modified diffusion equation** using an additional scalar field, $\gamma$, initialised in the droplet phase.
- **Quantitative mixing** is evaluated using a **custom mixing metric**, automatically written to a ```postProcessing/``` directory during simulation.

All solvers and utilities are based on the open-source finite volume framework [OpenFOAM 9](https://openfoam.org/), which must be installed prior to use.

## ⚙️ Problem Parameters
The droplet impact and mixing process can be characterized by the following key physical parameters:

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
Diffusivity                           | $D$           | system/controlDict
Sesile Droplet Volume                 | $V_s$         | system/setFieldsDict
Sessile Droplet Initial Contact Angle | $\theta_s$    | system/setFieldsDict

**Explain how the sessile droplet geometry is set through a psuedo radius etc.**

## 🧩 Installation
This code is compatible with OpenFOAM 9, available for download [here](https://openfoam.org/version/9/)

The repository mirrors the folder structure of a standard OpenFOAM installation, allowing it to be copied directly into your user directory (e.g. ```$HOME/OpenFOAM/user-9/```).

All custom solvers and utilities must be compiled using OpenFOAM’s build system:
```
wmake
```
Alternatively, after cloning the repository:
```
git clone https://github.com/yatin-darbar/Droplet_Mixing.git
cd Droplet_Mixing
wmake -all
```

### 🧪 Testing the Installation
To verify successful compilation, navigate to an example case (**e.g. examples/case1)** and run:
```
blockMesh
setFields
diffusiveInterFoam
```
If the solver runs without warnings in the terminal output, the installation is complete.

## 🧱 Mesh Generation
Mesh generation and field initialization are handled using the provided helper script:
```refineSetFields.sh```

**Note:** Ensure the script is executable (``` chmod +x refineSetFields.sh ```)

This script:
- Generates the initial simulation mesh
- Sets the droplet geometry and field distributions

It is good practice to maintain a ```0.orig``` directory that stores the pre-processed initial conditions (before ```refineSetFields.sh``` is run).
  
 Simulation can be started with the command:
```
diffusiveInterFoam
```


##  🧊 Static Diffusion
In cases where droplets have reached equilibrium and only diffusion needs to be simulated, the momentum equations can be deactivated to reduce computational cost.

Use the ```staticDiffusiveInterFoam``` solver to solve only the phase-restricted diffusion equation without the Navier–Stokes terms.

If droplets are stationary, the ```staticSetFieldsDict``` file can be used in place of the standard setFieldsDict to initialize a new simulation with:
- Unaltered droplet geometry, and
- A zero velocity field.

This allows the diffusion of the scalar field to evolve independently of flow dynamics.

## 📁 Repository Structure

Folder           | 	Description
--------         |---------------
```solvers/```   |	Source code for custom solvers
```utilities/``` |	Source code for the Kistler model and custom function objects
```examples/```	 | Example case files
```data/```	     | Data used in associated publications
```figures/```   |	Plots and figures from the study


## 🧾 Citation
If you use or build upon this solver, please cite:
```
@software{YatinDarbar_diffusiveInterFoam_2025,
  author = {Yatin Darbar},
  title = {Droplet Mixing -- GitHub Repository},
  year = {2025},
  url = {https://github.com/yatin-darbar/Droplet_Mixing}
}
```
(This will be updated with the publication citation once available.)


