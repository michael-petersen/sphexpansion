# sphexpansion
## libraries to run potential expansions 

This project provides a standalone interface to evaluate basis function expansion forces for a flexible number of components.
The utilities are headers only, so no compilation is required beyond your application.

Each component is described by four objects:
1. Sturm-Louiville solution basis functions
2. The spherically-symmetric monopole model
3. The weights for the basis through time
4. (optional) The location of the basis centre through time

The libraries to realise the potentials are in include/. Some working examples are located in the examples/ folder.

There are also working utilities:
1. A leapfrog integrator for test integrations
2. A (draft) method to accumulate expansion coefficients.
