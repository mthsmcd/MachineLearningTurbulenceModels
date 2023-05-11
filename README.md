# MachineLearningTurbulenceModels
OpenFOAM's Turbulence Models to be used with Machine Learning predictions.

These models were developed and tested in OpenFOAM-4.x and OpenFOAM-7

Models are used to correct RANS turbulent flows by using Machine Learning predicted quantities or direct injection of DNS
or other high-fidelity simulations'data (e.g. LES).

Models' source terms:
- **RStress** - Reynolds Stress Tensor ***R***
  - Directly injects the deviatoric part of ***R*** into the momentum balance
- **nutRStress** - Perpendicular-to S (mean strain-rate tensor) Reynolds Stress ***Rperp*** and optimal turbulent viscosity ***nut***
  - Based on the papers by Wu et al. (2018) and Brener et al. (2021), referenced at the end of this file.
  - Directly injects the deviatoric part of ***Rperp*** into the momentum balance, along with the turbulent viscosity nut.
  - The turbulent viscosity ***nut*** can be included in the diffusive term of the momentum equation implicitly, as proposed by Wu et al. (2018), or explicitly, as in Brener et al. (2021)
  - To select between implicit or explicit variations, a constant ***implicitFactor*** needs to be defined as a model coefficient in the `turbulenceProperties` dictionary. The constant needs to be assigned a value of `0.0` (explicit) or `1.0` (implicit).
  - If not defined, the simulation runs with the default value of `1.0` 
- **tForce** - Modified Reynolds Force Vector ***t***
  - Based on the work by Cruz et al. (2019), referenced at the end of this file.
  - Directly injects the vector ***t*** into the momentum balance
- **nutTForce** - nonlinear part of the Modified Reynolds Force Vector ***tStar*** and an optimal turbulent viscosity ***nut***
  - Based on the papers Brener et al. (2021), Brener et al (2022) and  Cruz et al (2019), referenced at the end of this file.
  - Directly injects the vector ***tStar*** into the momentum balance along the turbulent viscosity nut
  - The scalar ***nut*** is included within the diffusive term of the discretized mean momentum balance solved to compute the velocity field U
  - Analogous to the `nutRStress` model, a constant ***implicitFactor*** defines if the diffusive term containing ***nut*** is calculated implicitly or explicitly.
  - Default value is also `1.0` (implicit)

Models were constructed using OF's *ShihQuadraticKE* turbulence model.

To include the library in your OF installation use the command:
1) Pull the repository, preferably into your $WM_PROJECT_USER_DIR
2) Go to the directory where you copied the repository's content
3) Use the command `wmake libso`
4) To use the models it's necessary to include the line below into your simulation's controlDict:
  `libs ("libmyMachineLearningRASModels.so");`
5) Change the turbulence model in `constant/turbulenceProperties` into one of the 4 models of this library.
  
# References
- Brener, B. P., Cruz, M. A., Macedo, M. S. S. and Thompson, R. L. "An Invariant and Highly–Accurate Strategy for Data-Driven Turbulence Modelling." *SSRN Electronic Journal* (2022) http://dx.doi.org/10.2139/ssrn.4073177

- Brener, B. P., Cruz, M. A., Thompson, R. L., & Anjos, R. P. "Conditioning and accurate solutions of Reynolds average Navier–Stokes equations with data-driven turbulence closures."   *Journal of Fluid Mechanics*, 915, A110 (2021). https://doi.org/doi:10.1017/jfm.2021.148

- Cruz, M. A., Thompson, R. L., Sampaio, L. E., & Bacchi, R. D. "The use of the Reynolds force vector in a physics informed machine learning approach for predictive turbulence modeling." *Computers & Fluids* 192 (2019): 104258. https://doi.org/10.1016/j.compfluid.2019.104258

- Wu, J.L., Xiao, H., and Paterson, E. "Physics-informed machine learning approach for augmenting turbulence models: A comprehensive framework." *Physical Review Fluids* 3.7 (2018): 074602. https://doi.org/10.1103/PhysRevFluids.3.074602
