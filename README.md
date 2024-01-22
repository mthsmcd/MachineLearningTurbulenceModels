# MachineLearningTurbulenceModels
OpenFOAM's (OF) Turbulence Models to be used with Machine Learning predictions.

**These are the models used in our paper *"A highly accurate strategy for data-driven turbulence modeling"* by Bernardo P. Brener, Matheus A. Cruz, Matheus S. S. Macedo and Roney L. Thompson., published at *Computational and Applied Mathematics*** 

The paper can be accessed at https://doi.org/10.1007/s40314-023-02547-9, you can find information on how to cite it following this link.
It is fully available for free at https://rdcu.be/dwtcd

(It was previously cited by other works in its *preprint* version, which is still accessible at http://dx.doi.org/10.2139/ssrn.4073177)

**The models are used to correct RANS simulations by using quantities predicted by Machine Learning techniques or by the direct injection of high-fidelity fields (e.g. DNS, LES).**

The corrections are driven by source-terms injected into the mean momentum equation.

Implementation and tests were done in OpenFOAM-4.x, OpenFOAM-7 and OpenFOAM-2306.

*The OpenFOAM foundation versions (openfoam.org) have renamed and moved header files used to compile this library from version 8 onwards. In these versions, compilation won't succeed, unless the code is adapted. For this reason I advise anyone interested in using this library to prefer the ESI versions (openfoam.com)*

## Folders in the repository

### The folder `of-turbulence-models` contains the OpeFOAM implementation of the data-driven turbulence models.

To compile and include the library in your OF installation do the following:
1) Pull the repository, preferably into your $WM_PROJECT_USER_DIR
2) Navigate to the repository's directory `of-turbulence-models`
3) Use the command `wmake libso`
4) It's necessary to include the line
   `libs ("libMachineLearningTurbulenceModels.so");` into your simulation's controlDict in order to use the models
5) Change the turbulence model in `constant/turbulenceProperties` into one of the 4 models of this library.

### The folder `of-applications` contains the applications that calculate the source-terms of each model

To compile the applications, navigate to their directories and use the command `wmake`.
A shell script that compiles all of them is provided.

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

The DNS fields for the square-duct were provided by Pinelli et al. (2010) and post-processed by Fonseca et al. (2022).

The DNS fields for the periodic-hills were provided by Xiao et al. (2020)

## Models' in the repository and their source terms
- **RST** - Reynolds stress tensor ***R***
  - Directly injects the deviatoric part of ***R*** into the momentum balance
- **evRST** - Perpendicular-to S (mean strain-rate tensor) Reynolds Stress ***Rperp*** and optimal eddy-viscosity ***nut***
  - Based on the paper by Wu et al. (2018), referenced at the end of this file.
  - Directly injects the deviatoric part of ***Rperp*** into the momentum balance, along with the eddy-viscosity nut.
  - The turbulent viscosity ***nut*** can be included in the diffusive term of the momentum equation implicitly, as proposed by Wu et al. (2018), or explicitly, as in Brener et al. (2021)
  - Optionally, you can select between implicit or explicit variations, a constant ***implicitFactor*** can be defined as a model coefficient in the `turbulenceProperties` dictionary. If defined, the constant needs to be assigned a value of `0.0` (explicit) or `1.0` (implicit).
  - If not defined, the model assigns the default value of `1.0` 
- **RFV** - Modified Reynolds force vector ***t***
  - Based on the work by Cruz et al. (2019), referenced at the end of this file.
  - Directly injects the vector ***t*** into the momentum balance
- **evRFV** - nonlinear part of the modified Reynolds force vector ***tStar*** and an optimal eddy-viscosity ***nut***
  - Based on the papers Brener et al. (2021) and Brener et al (2022), referenced at the end of this file.
  - Directly injects the vector ***tStar*** into the momentum balance along with the eddy-viscosity nut.
  - The scalar ***nut*** is included within the diffusive term of the discretized mean momentum balance solved to compute the velocity field U
  - Analogous to the `RST-EV` model, a constant ***implicitFactor*** defines if the diffusive term containing ***nut*** is calculated implicitly or explicitly.
  - Default value is also `1.0` (implicit)

Models were constructed using OF's *ShihQuadraticKE* turbulence model.

Inside the `data` folder there is a shell script that will calculate and organize the source-terms in the simulations folders.


## References

**Models**

- Brener, B. P., Cruz, M. A., Macedo, M. S. S. and Thompson, R. L. "A highly accurate strategy for data-driven turbulence modeling." *Computational and Applied Mathematics*, 43, 59 (2024). https://doi.org/10.1007/s40314-023-02547-9

- Brener, B. P., Cruz, M. A., Thompson, R. L. and Anjos, R. P. "Conditioning and accurate solutions of Reynolds average Navier–Stokes equations with data-driven turbulence closures." *Journal of Fluid Mechanics*, 915, A110 (2021). https://doi.org/doi:10.1017/jfm.2021.148

- Cruz, M. A., Thompson, R. L., Sampaio, L. E. and Bacchi, R. D. "The use of the Reynolds force vector in a physics informed machine learning approach for predictive turbulence modeling." *Computers & Fluids*, 192 (2019): 104258. https://doi.org/10.1016/j.compfluid.2019.104258

- Wu, J.L., Xiao, H., and Paterson, E. "Physics-informed machine learning approach for augmenting turbulence models: A comprehensive framework." *Physical Review Fluids*, 3.7 (2018): 074602. https://doi.org/10.1103/PhysRevFluids.3.074602


**Databases**

- Pinelli, A., Uhlmann, M., Sekimoto, A. and Kawahara, G. "Reynolds number dependence of mean flow structure in square duct turbulence." *Journal of Fluid Mechanics*, 644, 107-122 (2010). https://doi.org/10.1017/S0022112009992242

- Xiao, H., Wu, J.-L., Laizet, S., Duan, L. "Flows over periodic hills of parameterized geometries: A dataset for data-driven turbulence modeling from direct simulations." *Computers and Fluids*, 200(104431), 1–26 (2020). https://doi.org/10.1016/j.compfluid.2020.104431

- Fonseca, E.F., Rangel, V.B., Brener, B.P., Cruz, M. A. and Thompson, R. L. "Pre-processing DNS data to improve statistical convergence and accuracy of mean velocity fields in invariant data-driven turbulence models." *Theoretical and Computational Fluid Dynamics*, 36, 435–463 (2022). https://doi.org/10.1007/s00162-022-00603-4
