# UniverseModels
Models to describe the evolution of the universe using EFT, dark energy and modified gravity theories. This repo contains four subfolders:

* ### UniverseModels
  This folder contains the necessary algorithms to model the background energy for any flat universe. This has been written in a Python notebook, therefore in order to open it you will need to have iPython installed (see http://ipython.org/install.html for more detail). However a python script of the code has also been included.
  * This only contains one script: WokingAndNomalised.ipynb, which when run plots the evolving scale factor with time and computational error.

* ### MinimizingEFT 
  This folder contains the algorithm for creating the χ2 models and mini- mizing these for multiple models. This contains four files:
  * The Main file: this run the entire script when compiled. This also contains the initial conditions for the models and defines which parameters are being fitted.
  * The Standard file: this evolves the ΛCDM case which is required in calculating χ2.
  * The ChiSq file: calculates the χ2 values for a model, given a specific scalar field al- teration. The scalar field being altered must be defined in this class, in the function definition.
  * The ChiSq4params file: calculates the χ2 values when all three scalar fields are being altered.

* ### Perturbed EFT 
  Allows individual parameters to be altered, as well as calculating the gravitational potential perturbations. This contains three files:
  * The Main file: this run the entire script when compiled. This also contains the initial conditions for the models and defines the scalar field parameters.
  * The Standard file: this evolves the ΛCDM case which is required in calculating χ2.
  * The ScalarPerturbations file: calculates the φ, ψ and δρ at each time-step.

* ### CurrentPerturbations
  Contains the code for calculating the perturbations of the EFT scalar fields. This also contains three files:
  * The Main file: this run the entire script when compiled. This also contains the initial conditions for the models and defines the scalar field parameters.
  * The Standard file: this evolves the ΛCDM case which is required in calculating χ2.
  * The ScalarPerturbations file: calculates the φ, ψ, π and δρ at each time-step.
