# UniverseModels
Models to describe the evolution of the universe using EFT, dark energy and modified gravity theories. This repo contains four subfolders:

* ### BackgroundEnergy
  This has been written in a Python notebook, therefore in order to open it you will need to have iPython installed (see http://ipython.org/install.html for more detail). Firstly, the background energy for a flat universe was modelled for a matter only, vacuum energy only and ΛCDM case (Einstein-de Sitter, deThis folder contains the necessary algorithms to model the background energy for any flat universe.  Sitter and concordance model respectively). This was done by evolving the Friedmann equation discretely as described. These theoretical results were compared to the analytical expansion. The errors were then found by calculating the difference between these methods at each time-step. These computations were conducted using a time-step of ∆t = 0.001Gyr and a system run for 20,000 iterations.

  * This only contains one script: WokingAndNomalised.ipynb, which when run plots the evolving scale factor with time and computational error.

* ### MinimizingEFT 
  This folder contains the algorithm for creating the χ2 models and mini- mizing these for multiple models. This contains four files:
  * The Main file: this run the entire script when compiled. This also contains the initial conditions for the models and defines which parameters are being fitted.
  * The Standard file: this evolves the ΛCDM case which is required in calculating χ<sup>2</sup>.
  * The ChiSq file: calculates the χ<sup>2</sup> values for a model, given a specific scalar field alteration. The scalar field being altered must be defined in this class, in the function definition.
  * The ChiSq4params file: calculates the χ<sup>2</sup> values when all three scalar fields are being altered.

* ### PerturbedEFT 
  Allows individual parameters to be altered, as well as calculating the gravitational potential perturbations. This contains three files:
  * The Main file: this run the entire script when compiled. This also contains the initial conditions for the models and defines the scalar field parameters.
  * The Standard file: this evolves the ΛCDM case which is required in calculating χ<sup>2</sup>.
  * The ScalarPerturbations file: calculates the φ, ψ and δρ at each time-step.

* ### CurrentPerturbations
  Contains the code for calculating the perturbations of the EFT scalar fields. This also contains three files:
  * The Main file: this run the entire script when compiled. This also contains the initial conditions for the models and defines the scalar field parameters.
  * The Standard file: this evolves the ΛCDM case which is required in calculating χ<sup>2</sup>.
  * The ScalarPerturbations file: calculates the φ, ψ, π and δρ at each time-step.
