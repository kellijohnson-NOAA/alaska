# TODO

1. gstat: a spatial r package that will facilitate mapping of spatial predictions
using a inverse linear distance kriging algorithm.
2. Plot the mean log abundance for each stock behind the plot for the overall mean log abundance,
I am going to have to change what is placed in the report to get at the predicted
abundance for each station in a given year.
3. Use VAST
4. @Mueter2005 also states why years prior to 1993 should be removed
because of a lack of standardized gear, and they suggest that
only trawls with a performance of satisfactory should be included, but I am
not sure how to find out which trawls are and are not satisfactory.
5. Rerun the analysis with the SE Alaska no-trawl zone kept in
6. Read: Lindgren F, Rue H, Lindström J. 2011.
An explicit link between Gaussian fields and Gaussian Markov random fields:
the stochastic partial differential equation approach [with discussion]. J. R. Stat. Soc. B. 73(4):423–98.
7. Spatial clusters were attempted with a multidirectional optimum ecotope-based algorithm,
which uses nearest neighbors and a dependent variable to determine clusters.
The following code will run the analysis, but it is dependent on how far
you assume a neighbor can be from another point, and was thus not investigated any further.
# library(AMOEBA)
# res <- AMOEBA(omega, dnearneigh(coordinates(info), 0, 0.10), 1, 1)
8. Find a correction factor to downweight the kmeans clustering, where the correction factor should be based on the estimated variance and the difference between the true differences in bettering the model with an increase in cluster number when there really are two true clusters and when there are no true clusters. This will lead to some clusters potentially being missed when they do exist, but maybe a decrease in the number of clusters being identified when there are none. Suggested by JTT.
9. Table order
10. Figure order
11. NOAA JISAO contribution number.

# Notes

## 20140426 Kerim Aydin
* Many time-series analyzes are available for the Bering Sea, but an analysis of the Gulf of Alaska are fewer. Therefore, the study will be a good counter-point to available analyses.
* AIs have a strong longitudinal gradient
* AIs house the western stock of stellar sea lions, which are endangered, while the eastern stock is no longer endangered.
* It would be good to look at community composition trends and see if and how they have changed over time. Furthermore, if there are trends between groups as well as trends within groups. In the AIs, cod is more of a keystone predator than in other areas, making for a situation with greater effects for the same level of a perturbation.
* Quota management for cod in the AIs is highly contentious

## References
* Ivonne Ortiz - dissertation looked at the food web in 2 degree blocks
* Goundie, E. T., David, A. S., and Trites, A. W. 2015. Low prey abundance leads to less efficient foraging behavior in Steller sea lions. Journal of Experimental Biology and Ecology, 470:70-77.
* Hui, T. C. Y., Gryba, R. Gregr, J., and Trites, A. W. 2015. Assessment of competition between fisheries and Steller sea lions in Alaska based on estimated prey biomass, fisheries removals and predator foraging behaviour. PLOS ONE 10(5), e0123786.
* Jansen, T., Kristensen, K., Kainge, P., Durholtz, D., Strømme, T., Thygesen, U. H., Wilhelm, M. R., Kathena, J., Fairweather, T. P., Paulus, S., Degel, H., Lipinski, M. R., and Beyer, J. E. 2016. Migration, distribution and population (stock) structure of shallow-water hake (*Merluccius capensis*) in the Benguela Current Large Marine Ecosystem inferred using a geostatistical population model. Fisheries Research, 179: 156-167.
* Rushing, C. S., Ryder, T. B., Scarpignato, A. L., Saracco, J. F., and Marra, P. P. 2015. Using demographic attributes from long-term monitoring data to delineate natural populations structure. Journal of Applied Ecology
