# sphexpansion (MWLMC)
### Milky Way - Large Magellanic Cloud simulation from Lilleengen et al. (2022)

This branch of sphexpansion provides interfaces a model Milky Way and Large Magellanic Cloud interaction.

------------------
(user interface to be described here)

------------------

### Installation notes

sphexpansion now uses `cmake` as the primary installation tool (to build the example applications). To use, create a directory `build`, then in that directory, `cmake ..; make` to create working examples in `build/examples`. For compiling with `pybind11` capability, you may want to use `cmake -DPYTHON_EXECUTABLE=$(python3 -c "import sys; print(sys.executable)") ..; make`.

You may need to install `cmake`. Using homebrew, try `brew install cmake`. If using macports, the corresponding command is `sudo port install cmake`.

You may also need to install `eigen` (v3), a header-only C++ library. As a header only library, it is sufficient to download the tarball, unpack, and copy the `Eigen/` directory to `/usr/local/include` (which is already part of the search path for clang). You may also use `brew install eigen` if using a Mac with homebrew. If using macports, the command is `sudo port install eigen3`.


-----------------------------

### Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
