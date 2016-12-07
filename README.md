# alaska
Using spatially-explicit Gompertz population dynamics models to estimate the 
population structure of Pacific cod *Gadus macrocephalu* in the Aleutian Islands (AI) and
the Gulf fo Alaska (GOA).
The analysis includes a simulation experiment where a theoretical fish population with
either one or two subpopulations is generated and fit to both the Gompertz model and
a spatial algorithm to estimate subpopulation structure based on local productivity.

## Reference
Currently, the manuscript is mid stage and will eventually be submitted to 
the Journal of Applied Ecology.
Johnson *et al.* in prep. Using gradients in productivity to inform population
structure for a managed species: a case study using Pacific cod in Alaska.
Journal of Applied Ecology 00: 00-00.

## To do list
- [ ] generate a distance matrix with space and estimates of productivity
- [ ] create a function that will estimate the management unit using 
clusters, such as `cluster::pam` or `pamk`
- [ ] resample omega using the standard deviation and a mean of zero
- [ ] look at how changing the mesh boundary changes estimates, 
according to the documentation the mesh boundary needs to be at least equal
to the range to avoid bias
- [ ] look at using something other than a normal distribution for the 
continuous estimate of omega in the spatial algorithm
- [x] email JTT about the lack of bias removal because of changes
to the recursive `.cpp` file.
