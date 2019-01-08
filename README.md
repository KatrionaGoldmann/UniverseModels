# UniverseModels
Models to describe the evolution of the universe using EFT, dark energy and modified gravity theories.

This repo contains four subfolders, containing generally one Python notebook (so in order to open it you will need to have iPython installed - see http://ipython.org/install.html for more detail) and necessary python scripts:

----
* ### 1. BackgroundEnergy

  The script in this directory models the background energy for a flat universe looking at

    i) a matter only (Einstein-de Sitter),

    ii) vacuum energy only (deSitter),

    iii) and ΛCDM case (and concordance model).

  This is achieved by evolving the Friedmann equation discretely. The theoretical results are compared to the analytical/computational expansion. The errors are then found by calculating the difference between these methods at each time-step. These computations are, by default, conducted using a time-step of ∆t = 0.001Gyr and a system run for 20,000 iterations.

#### Scripts:
  * BackgroundEnergy.ipynb, which when run plots the evolving scale factor with time and computational error.

----
* ### 2. PerturbedEFT

The analysis finds the overdensities evolving in EdS, dS and MG models under the Newtonian gauge. This is achieved considering the case without any anisotropic pressure, therefore φ = ψ. The plots below show the scale factor a(t), gravitational potential φ, and overdensity δ evolving over time. These were evolved over time increments using 50,000 iterations.

#### Scripts:
  * PerturbedEFT.ipynb: this run the entire script when compiled. This contains the initial conditions for the models and defines the scalar field parameters.
----
* ### 3. MinimizingEFT
  This folder contains the algorithm for creating the χ<sup>2</sup> models and minimising these for multiple models.

#### Scripts:
  * MinimizingEFT.ipynb: this run the entire script when compiled. This also contains the initial conditions for the models and defines which parameters are being fitted.
  * The EFT file: calculates the φ, ψ, π and δρ at each time-step.
  * The ChiSq file: calculates the χ<sup>2</sup> values for a model, given a specific scalar field alteration. The scalar field being altered must be defined in this class, in the function definition.
  * The ChiSq4params file: calculates the χ<sup>2</sup> values when all three scalar fields are being altered.

----

* ### CurrentPerturbations
  Contains the code for calculating the perturbations of the EFT scalar fields.

#### Scripts:
  * CurrentPerturbations.ipynb: this run the entire script when compiled. This also contains the initial conditions for the models and defines the scalar field parameters.
  * The ChiSq4params file: this evolves the ΛCDM case which is required in calculating χ<sup>2</sup>.
  * The EFT file: calculates the φ, ψ, π and δρ at each time-step.
