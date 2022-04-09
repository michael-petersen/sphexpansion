# TODO
## List of tasks to fill out this branch

0. Strip out
  - have velocity flag
  - duplicated double definitions?
  - all flags will be hardcoded
1. Convert to eigen (remove boost)
  - Need 3d strategy for coefficients, which should be vector of MatrixXd, probably.
2. Hard code loading up MW and LMC spheres (load on initialisation)
  - Add a default static MW, LMC coefficient file
3. Hard code loading up disc
4. Remove examples (replace!)
5. Strip comments
6. Add -fopenmp support, if available
7. Overload several transformation functions with eigen operations
