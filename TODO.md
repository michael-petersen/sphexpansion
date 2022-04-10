# TODO
## List of tasks to fill out this branch

0. Strip out
  - havevelocity flag
  - duplicated double definitions?
  - all flags will be hardcoded
1. Convert to eigen (remove boost)
  - Need 3d strategy for coefficients, which should be vector of MatrixXd, probably.
2. Hard code loading up MW and LMC spheres (load on initialisation)
  - Add a default static MW, LMC coefficient file
3. Hard code loading up disc
4. Replace examples
  - Check rotation curves to start
5. Strip comments
6. Add -fopenmp support, if available (modern computers only)
7. Overload several transformation functions with eigen operations
