# MachineLearningTurbulenceModels
OpenFOAM implementation of turbulence models driven by Machine Learning predictions.

**These models were used in our papers:**

1) ***"A highly accurate strategy for data-driven turbulence modeling"* by Bernardo P. Brener, Matheus A. Cruz, Matheus S. S. Macedo and Roney L. Thompson., published at *Computational and Applied Mathematics* in January 2024** 

    Available at https://doi.org/10.1007/s40314-023-02547-9 , you can read it for free [on this link](https://rdcu.be/dwtcd). 

    (It was previously cited by other works in its *preprint* version, still accessible [here](http://dx.doi.org/10.2139/ssrn.4073177))


2) ***"A data-driven turbulence modeling for the Reynolds stress tensor transport equation"* by Matheus S. S. Macedo, Matheus A. Cruz, Bernardo P. Brener and Roney L. Thompson, published at the *International Journal for Numerical Methods in Fluids* in March 2024**

    Available at https://doi.org/10.1002/fld.5284 , you can read it for free [on this link](https://onlinelibrary.wiley.com/share/author/JQGYVIFEU2NRDDMCZRJZ?target=10.1002/fld.5284)


Information on how to cite both papers are available on their respective links

**The models are used to correct RANS simulations by using quantities predicted by Machine Learning techniques. They can also be used by the directly using high-fidelity fields (e.g. DNS, LES).**

In the models of the 1st paper, the corrections are driven by source terms injected into the mean momentum equation, while the model of the 2nd paper injects its source term into a Reynolds stress model (RSM).


## Table of Contents

- [Compatibility](#compatibility)
- [Folders in the repository](#folders-in-the-repository)
  - [of-turbulence-models](#the-folder-of-applications-contains-the-applications-that-calculate-the-source-terms-of-each-model)
  - [of-applications](#the-folder-of-applications-contains-the-applications-that-calculate-the-source-terms-of-each-model)
  - [data](#the-folder-data-contains-openfoam-simulations)
- [Models in the repository and their source terms](#models-in-the-repository-and-their-source-terms)
- [Using the models](#using-the-models)
  - [Tutorials: using the models with the provided DNS data](#tutorials-using-the-models-with-the-provided-dns-data)
    - [Correcting the square-duct RANS case with Re = 3500 using the RST](#1---correcting-the-square-duct-rans-case-with-re--3500-using-the-rst)
    - [Correcting the periodic-hill RANS case with alpha = 1.0 using the EV-RFV](#2---correcting-the-periodic-hill-rans-case-with-alpha--10-using-the-ev-rfv)
    - [Correcting the square-duct RANS case with Re = 2200 using the Gamma-RST](#3---correcting-the-square-duct-rans-case-with-re--2200-using-the-gamma-rst)
- [References](#references)
  - [Models](#models)
  - [Databases](#databases)


## Compatibility

The models were implemented using the OF *ShihQuadraticKE* model as a base.
Implementation and tests were done in OpenFOAM-4.x, OpenFOAM-7, OpenFOAM-v2306 and OpenFOAM-v2312.

*The OpenFOAM foundation versions (openfoam.org) have renamed and moved header files used to compile this library from version 8 onwards. In these versions, compilation won't succeed, unless the code is adapted. For this reason I advise anyone interested in using this library to prefer the ESI versions (openfoam.com)*

## Folders in the repository

### The folder `of-turbulence-models` contains the OpeFOAM implementation of the data-driven turbulence models.

### The folder `of-applications` contains the applications that calculate the source terms of each model

To compile and use the applications and the libraries, you need to do the following:
1) Clone the repository, preferably into your $WM_PROJECT_USER_DIR
2) Navigate to the repository's directory and execute the scripts `./Allwclean` and `./Allwmake` 
4) After compilation, it is necessary to include the line
   `libs ("libMachineLearningTurbulenceModels.so");` into your simulation's controlDict in order to use the turbulence models
5) Change the turbulence model in `constant/turbulenceProperties` into one of the 5 models of this library.

### The folder `data` contains OpenFOAM simulations

The square-duct (SD) and periodic-hills (PH) simulations used in our paper are provided. 

The SD folder contains the simulations for Reynolds numbers of 2200, 2400, 2600, 2900, 3200, 3500.
The PH folder contains the simulations for the slopes of 0.5, 0.8, 1.0, 1.2, 1.5.
Each subdirectory in the SD or PH folders contains the subdirectory `0` with the following `k-epsilon` fields:
- *Urans* - velocity
- *Rrans* - Reynolds stress
- *p* - pressure
- *S* - mean strain-rate tensor
- *k* - turbulent kinetic energy
- *epsilon* - turbulent dissipation
- *nut* - eddy-viscosity

And the following DNS fields:
- *Udns* - velocity
- *Rdns* - Reynolds stress

The DNS fields for the square-duct were provided by [Pinelli et al. (2010)](#references) and post-processed by [Fonseca et al. (2022)](#references).

The DNS fields for the periodic-hills were provided by [Xiao et al. (2020)](#references)

## Models' in the repository and their source terms
- **RST** - Reynolds stress tensor ***R***
  - Directly injects the deviatoric part of ***R*** into the momentum balance
- **evRST** - Perpendicular-to S (mean strain-rate tensor) Reynolds Stress ***Rperp*** and optimal eddy-viscosity ***nut***
  - Based on the paper by [Wu et al. (2018)](#references)
  - Directly injects the deviatoric part of ***Rperp*** into the momentum balance, along with the eddy-viscosity nut
  - The turbulent viscosity ***nut*** can be included in the diffusive term of the momentum equation implicitly, as proposed by [Wu et al. (2018)](#references), or explicitly, as in [Brener et al. (2021)](#references)
  - You can select between implicit or explicit variations, a constant ***implicitFactor*** can be defined as a model coefficient in the `turbulenceProperties` dictionary. If defined, the constant needs to be assigned a value of `0.0` (explicit) or `1.0` (implicit).
  - If not defined, the model assigns the default value of `1.0` 
- **RFV** - Modified Reynolds force vector ***t***
  - Based on the work by [Cruz et al. (2019)](#references)
  - Directly injects the vector ***t*** into the momentum balance
- **evRFV** - nonlinear part of the modified Reynolds force vector ***tStar*** and an optimal eddy-viscosity ***nut***
  - Based on the papers by [Brener et al. (2021)](#references) and [Brener et al (2024)](#references)
  - Directly injects the vector ***tStar*** into the momentum balance along with the eddy-viscosity nut.
  - The scalar ***nut*** is included within the diffusive term of the discretized mean momentum balance solved to compute the velocity field U
  - Analogous to the `RST-EV` model, a constant ***implicitFactor*** defines if the diffusive term containing ***nut*** is calculated implicitly or explicitly.
  - Default value is also `1.0` (implicit)
- **gammaRST** - Symmetric source term tensor ***Gamma***
  - Based on the paper by [Macedo et al. (2024)](#references)
  - Injects the source term ***Gamma*** into the data-driven Reynolds stress model.
  - At each iteration, the RST equation is solved for ***R***.
  - The deviatoric part of the calculated ***R*** is injected into the momentum balance. 
  - The process is repeated iteratively until numerical convergence.

Inside the `data` folder there is a shell script that will calculate and organize the source terms in the simulations folders.

## Using the models

You can use the DNS fields contained in the `data` folder simulations to calculate the source terms of 
each turbulence model and employ them as the target of your ML.
The RANS data should be used to calculate the inputs of the ML technique you have selected.

The feature selection of the RANS inputs is itself a topic of discussion and, therefore, varies from 
work to work.
For example, almost each of the [referenced works](#references) uses a different set of inputs.
Because the focus of this repository is the OpenFOAM implementation of the turbulence models, we have
decided to only include the baseline RANS fields from where we extracted the features used as 
our inputs.

There are multiple manners of training ML techniques and obtaining the source terms of each of the 
5 models.
In our works we have used neural networks and random forests, they were both built using
common Python libraries such as: `Keras`, `TensorFlow` and `Scikit-learn`. 
For more details on the implementation and training of the techniques 
please check our papers and others in the [reference section](#references) of this README.

If your main interest is building ML architectures and predicting the fields, you are referred to the 
work of [McConkey et al (2021)](https://www.nature.com/articles/s41597-021-01034-2), where the authors
have curated an extensive database readily usable in ML algorithms, especially in Python.
The database is composed by DNS, LES and RANS data, and available at the
Kaggle platform [on this link](https://www.kaggle.com/datasets/ryleymcconkey/ml-turbulence-dataset). 
The data is presented both as OpenFOAM cases and as `Numpy` arrays, multiple derived RANS feature 
fields usable as ML inputs are also available.

### Tutorials: using the models with the provided DNS data

You can also correct the RANS cases of this repository by running the models with the DNS fields 
included here.
This serves as an upper performance estimate for each of the 5 methodologies, since the best possible
scenario for ML is to exactly predict its targets.

Instructions for running three of these evaluations are provided below.
**They must be followed after the models and applications of the repository have been compiled in your OpenFOAM installation!**

The first two tutorials with the `RST` and `evRFV` models should converge quickly, the last tutorial, using the `gammaRST`
model, takes longer since 6 additional PDEs for the Reynolds stress must be solved.

(To simplify some of the steps in the tutorials you may run the provided 
shell script in the `data` folder.)

#### 1 - Correcting the square-duct RANS case with Re = 3500 using the RST

The model solves for `U` and `p` and requires the source term `R`.

0) Navigate to the case folder: `cd data/square-duct/3500/`
1) Use the RANS velocity field `Urans` as the initial condition for `U`: `cp -r 0/Urans 0/U`
2) Use the DNS Reynolds stress tensor `Rdns` as the model's source term: `cp -r 0/Rdns 0/R`
3) Specify `RST` as the turbulence model in `constant/turbulenceProperties`: `foamDictionary constant/turbulenceProperties -entry RAS.RASModel -set RST`
4) Include the repository's library in `system/cotrolDict`: `foamDictionary system/controlDict -entry libs -add "("'"libMachineLearningTurbulenceModels.so"'")"`

    **You only have to do this if you are not using any of the `controlDict` included in `data` folder**

5) Run the simulation: `simpleFoam` (you may use any other steady-state solver you prefer)


#### 2 - Correcting the periodic-hill RANS case with alpha = 1.0 using the EV-RFV

The model solves for `U` and `p` and requires the source terms `nut` (optimal) and `tStar`.

0) Navigate to the case folder: `cd data/periodic-hills/1p0/`
1) Calculate the DNS eddy-viscosity `nut` and nonlinear part of the Reynolds force vector `tStar`

    1) Rename the baseline RANS eddy-viscosity to avoid overwriting it: `mv -v 0/nut 0/nutRans`
    2) Set `Udns` and `Rdns` as `U` and `R`
       1) `cp -r 0/Udns 0/U`
       2) `cp -r 0/Rdns 0/R`
    3) Calculate the DNS `nut`: `calculateNut`
    4) Calculate the DNS `tStar`: `calculateRFVperp`
    5) Reorganize the `0` folder: `rm -r 0/U 0/R`

2) Use the RANS velocity field `Urans` as the initial condition for `U`: `cp -r 0/Urans 0/U`
3) Specify `evRFV` as the turbulence model in `constant/turbulenceProperties`: `foamDictionary constant/turbulenceProperties -entry RAS.RASModel -set evRFV`
4) Include the repository's library in `system/cotrolDict`: `foamDictionary system/controlDict -entry libs -add "("'"libMachineLearningTurbulenceModels.so"'")"`

   **You only have to do this if you are not using any of the `controlDict` included in `data` folder**

5) Run the simulation: `simpleFoam` (you may use any other steady-state solver you prefer)


#### 3 - Correcting the square-duct RANS case with Re = 2200 using the Gamma-RST

The model solves for `U`, `p` and `R` and requires the source terms `nut` (RANS) and `Gamma`.

0) Navigate to the case folder: `cd data/square-duct/2200/`
1) Calculate the DNS source term `Gamma`

    1) Set `Udns` and `Rdns` as `U` and `R`
        1) `cp -r 0/Udns 0/U`
        2) `cp -r 0/Rdns 0/R`
    3) Calculate the DNS `Gamma`: `calculateGamma`
    4) Reorganize the `0` folder: `rm -r 0/U 0/R`

2) Use the RANS velocity field `Urans` as the initial condition for `U`: `cp -r 0/Urans 0/U`
3) Use the RANS Reynolds stress field `Rrans` as the initial condition for `R`: `cp -r 0/Rrans 0/R`
4) Set the null boundary condition for `R`:
   1) `foamDictionary 0/R -entry boundaryField.fixedWalls.type -set fixedValue`
   2) `foamDictionary 0/R -entry boundaryField.fixedWalls.value -set "uniform (0 0 0 0 0 0)"`
5) *Optional - convergence occurs faster when the initial condition for `R` is null:* `foamDictionary 0/R -entry internalField -set "uniform (0 0 0 0 0 0)"`
6) Specify `gammaRST` as the turbulence model in `constant/turbulenceProperties`: `foamDictionary constant/turbulenceProperties -entry RAS.RASModel -set gammaRST`
7) Include the repository's library in `system/cotrolDict`: `foamDictionary system/controlDict -entry libs -add "("'"libMachineLearningTurbulenceModels.so"'")"`

   **You only have to do this if you are not using any of the `controlDict` included in `data` folder**

8) Run the simulation: `simpleFoam` (you may use any other steady-state solver you prefer)


## References

### Models

- Macedo, M. S. S., Cruz, M. A., Brener, B. P. and Thompson, R. L. "A data-driven turbulence modeling for the Reynolds stress tensor transport equation" *International Journal for Numerical Methods in Fluids* (2024). https://doi.org/10.1002/fld.5284

- Brener, B. P., Cruz, M. A., Macedo, M. S. S. and Thompson, R. L. "A highly accurate strategy for data-driven turbulence modeling." *Computational and Applied Mathematics*, 43, 59 (2024). https://doi.org/10.1007/s40314-023-02547-9

- Brener, B. P., Cruz, M. A., Thompson, R. L. and Anjos, R. P. "Conditioning and accurate solutions of Reynolds average Navier–Stokes equations with data-driven turbulence closures." *Journal of Fluid Mechanics*, 915, A110 (2021). https://doi.org/doi:10.1017/jfm.2021.148

- Cruz, M. A., Thompson, R. L., Sampaio, L. E. and Bacchi, R. D. "The use of the Reynolds force vector in a physics informed machine learning approach for predictive turbulence modeling." *Computers & Fluids*, 192 (2019): 104258. https://doi.org/10.1016/j.compfluid.2019.104258

- Wu, J.L., Xiao, H., and Paterson, E. "Physics-informed machine learning approach for augmenting turbulence models: A comprehensive framework." *Physical Review Fluids*, 3.7 (2018): 074602. https://doi.org/10.1103/PhysRevFluids.3.074602


### Databases

- Pinelli, A., Uhlmann, M., Sekimoto, A. and Kawahara, G. "Reynolds number dependence of mean flow structure in square duct turbulence." *Journal of Fluid Mechanics*, 644, 107-122 (2010). https://doi.org/10.1017/S0022112009992242

- Xiao, H., Wu, J.-L., Laizet, S., Duan, L. "Flows over periodic hills of parameterized geometries: A dataset for data-driven turbulence modeling from direct simulations." *Computers and Fluids*, 200(104431), 1–26 (2020). https://doi.org/10.1016/j.compfluid.2020.104431

- Fonseca, E.F., Rangel, V.B., Brener, B.P., Cruz, M. A. and Thompson, R. L. "Pre-processing DNS data to improve statistical convergence and accuracy of mean velocity fields in invariant data-driven turbulence models." *Theoretical and Computational Fluid Dynamics*, 36, 435–463 (2022). https://doi.org/10.1007/s00162-022-00603-4
