# mwlmc
### Milky Way - Large Magellanic Cloud simulation from Lilleengen et al. (2022)

This repository provides an interface to a simulation of the Milky Way and Large Magellanic Cloud (LMC) interaction.

------------------
### Documentation

The main functionalities are explained in the tutorials: 
 - Running [orbits](tutorial/orbits_tutorial.ipynb) in the simulation with options for different potential setups and 2 examples of dwarf galaxies;
 - Accessing the simulation [fields](tutorials/fields_tutorial.ipynb) (i.e. densities, potentials, and forces) for different potential setups.

**This simulation focuses on the effects of the interaction of the Milky Way and LMC on their haloes. Due to the evolving disc, computations in the inner Milky Way might not be reliable.** An example of a Sun-like orbit in the deforming potential vs a rigid Milky Way is given in the orbit tutorial.

------------------

### Installation notes
```
git clone --recursive https://github.com/sophialilleengen/mwlmc.git
pip install -e ./mwlmc 
```

You may need to install [`eigen` (v3), a header-only C++ library](https://eigen.tuxfamily.org/dox/GettingStarted.html). As a header only library, it is sufficient to download the tarball, unpack, and copy the `Eigen/` directory to `/usr/local/include` (which is already part of the search path for clang). You may also use `brew install eigen` if using a Mac with homebrew. If using macports, the command is `sudo port install eigen3`.

You may also need to install `cmake`. Using homebrew, try `brew install cmake`. If using macports, the corresponding command is `sudo port install cmake`.


-----------------------------
### Attributions
If you make use of this code, please cite the following papers
```
@ARTICLE{EXP,
       author = {{Petersen}, Michael S. and {Weinberg}, Martin D. and {Katz}, Neal},
        title = "{EXP: N-body integration using basis function expansions}",
      journal = {\mnras},
         year = 2022,
        month = mar,
       volume = {510},
       number = {4},
        pages = {6201-6217},
          doi = {10.1093/mnras/stab3639}
}
```

and consider citing this paper

```
@ARTICLE{Weinberg1999,
       author = {{Weinberg}, Martin D.},
        title = "{An Adaptive Algorithm for N-Body Field Expansions}",
      journal = {\aj},
         year = 1999,
        month = jan,
       volume = {117},
       number = {1},
        pages = {629-637},
          doi = {10.1086/300669}
}
```

-----------------------------

### License

Copyright 2022 Sophia Lilleengen and Michael Petersen.

mwlmc is free software made available under the MIT License. For details see the [LICENSE](./LICENSE) file.

-----------------------------

### Authors

Mike Petersen -  @michael-petersen - petersen@iap.fr \
Sophia Lilleengen - @sophialilleengen - s.lilleengen@surrey.ac.uk
