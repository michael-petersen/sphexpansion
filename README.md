# sphexpansion
### libraries to run potential expansions

This project provides a standalone interface to evaluate basis function expansion forces for a flexible number of components.
The main libraries are headers only, so no compilation is required beyond your application.
For Python extensions, compilation may be required.

Each component is described by four objects:
1. Sturm-Louiville solution basis functions
2. The spherically-symmetric monopole model
3. The weights for the basis through time
4. (optional) The location of the basis centre through time

The libraries to realise the potentials are in include/. Some working examples are located in the examples/ folder.

------------------

### Notes

sphexpansion now uses `cmake` as the primary installation tool (to build the example applications). To use, create a directory `build`, then in that directory, `cmake ..; make` to create working examples in `build/examples`.

You may have to install `eigen` (v3), a header-only C++ library. As a header only library, it is sufficient to download the tarball, unpack, and copy the `Eigen/` directory to `/usr/local/include` (which is already part of the search path for clang). You may also use `brew install eigen` if using a Mac with homebrew. If using macports, the command is `sudo port install eigen3`.

-----------------------------

### Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
