# 2D_DistanceSampling
 
In distance sampling surveys, the animals might avoid both the transects in the absence of observers, and the observers themselves. To correct for the effect of the behavioral responses of the animals to either the transects or the observers, we can estimate density and abundance using line transect survey data with both the forward and perpendicular distances to the observers (2D distance sampling - R LT2D package, Borchers and Cox 2015), not just the perpendicular distance. This analysis approach was also applied and recommended by Elenga et al. (2020).

Here, we rely on the functions from LT2D package (https://github.com/david-borchers/LT2D), as partly revised and applied in Elenga et al. (2020) (https://github.com/cbonenfant/duikers-abundance). With respect to the latter, we made additional minor changes to the code.

As an example, we apply the code to estimate the density of two species with different behavioral response to the transects/observers - impala and duikers.


## Density Surface Modelling


## Instructions

**Please also note that we assume that the folder in which this .Rmd file is stored includes the subfolders named *data* and *output*. The first one should contain the input data, while the second one must be created as an empty folder - output of the DS analyses will be saved there for (optional) subsequent analyses. The subfolder *functions* is also essential since it contains the customized functions used for the analyses. The whole structure can however be recreated effortless by forking and then cloning the GitHub repository on your local machine.**