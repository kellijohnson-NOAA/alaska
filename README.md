# alaska
Using spatially-explicit Gompertz population dynamics models to estimate the
population structure of Pacific cod *Gadus macrocephalu* in the Aleutian Islands (AI) and
the Gulf of Alaska (GOA).
The analysis includes a simulation experiment where a theoretical fish population with
either one or two subpopulations is generated and fit to both the Gompertz model and
a spatial algorithm to estimate subpopulation structure based on local productivity.

## Reference
Johnson *et al.* in prep. Using gradients in productivity to inform population
structure for a managed species: a case study using Pacific cod in Alaska.
Journal of Applied Ecology 00: 00-00.

## To do list
- [ ] Fit the subpopulation algorithm to the true value of omega
- [ ] Parameters have priors, or are started from a given value,
these starting values need to be investigated and how they change or do not
change parameter estimates. JTT used `rnorm` to start vector of random
effects, where others have used zero as a starting value
- [ ] resample omega using the standard deviation and a mean of zero
- [ ] look at how changing the mesh boundary changes estimates,
according to the documentation the mesh boundary needs to be at least equal
to the range to avoid bias
- [ ] email JTT in response to the bias in parameter estimates
- [ ] Change the .cpp file to use more of the INLA functions
- [ ] Estimate separate processes for the spatial effects rather than
sharing parameters such as `kappa` and `Q`
- [ ] Investigate bias in the estimates of spatial processes in the mesh
border versus main points inside the mesh
- [ ] Use variograms to estimate the spatial correlation in random effects
- [ ] Recreate the Lindgren *et al*. (2011) study to determine how parameters
are estimated within the .cpp.
- [ ] Would the subpopulation estimator work on process error, which is estimated
on a yearly basis, better than on spatial variation in productivity?
- [ ] Fix parameters at their true value and determine how parameter estimates change,
are there some parameters that are harder to estimate than others?
