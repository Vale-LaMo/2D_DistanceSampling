# 2D_DistanceSampling
 
In distance sampling surveys, the animals might avoid both the transects in the absence of observers, and the observers themselves. To correct for the effect of the behavioral responses of the animals to either the transects or the observers, we can estimate density and abundance using line transect survey data with both the forward and perpendicular distances to the observers (2D distance sampling - R LT2D package, Borchers and Cox 2017), not just the perpendicular distance. This analysis approach was also applied and recommended by Elenga et al. (2020).

Here, we rely on the functions from LT2D package (https://github.com/david-borchers/LT2D), as partly revised and applied in Elenga et al. (2020) (https://github.com/cbonenfant/duikers-abundance). With respect to the latter, we made additional minor changes to the code.

As an example, we apply the code to estimate the density of two species with different behavioral response to the transects/observers - impala and duikers.

## Density Surface Modelling

We follow up on the distance sampling analyses in two dimensions, and we use the corrected detection function, taking into account the behavioral response of the animals, to fit a density surface model (DSM, Hedley & Buckland 2004, Miller et al. 2013) to the data. Indeed, DSMs follow a two-stage approach: after accounting for detectability via distance sampling methods, they then model distribution via a generalized additive model.
While 2D distance sampling accounts for behavioral response, it does not overcome the problem of non-random transect placement, which is common to most field surveys, especially in rugged terrain. The DSM approach allows to correct for non-random transect placement, modelling the abundance (or density) of animals as a function of spatially explicit covariates, while also adjusting by the global detection probability. The expected abundance is related to the spatial covariates using a General Additive Model GAMs (Hastie & Tibshirani 1990). Possible predictors tipically include UTM coordinates, habitat, distance from water bodies, ecc.
The model is then used to produce a prediction grid of animal abundance within the study area. 


## Instructions

In order to replicate the analyses, or to apply them to your datasets, fork this repository and then clone it to your machine. This will re-create the working environment, with the main file *2D_DistanceSampling.Rproj* to be opened in RStudio.
The code for the analyses is stored in the *.Rmd* files *00_LT2D* and *01_dsm*.    
The input data are in the *data* folder.
The *functions* folder includes the scripts with the (customized) R functions.   
The *output* folder is meant to store the results of the analyses (mainly fitted models but also datasets produced in the first stage of the analysis and then used for DSMs).   

The list of packages needed for the analyses is reported at the beginning of each *.Rmd* file. Please make sure to install them before running the code. The LT2D package should be installed via devtools.


## References

Borchers DL, Cox MJ (2017) Distance sampling detection functions: 2D or not 2D? Biometrics, 73(2):593-602    
Hastie TJ, Tibshirani RJ (1990) Generalized Additive Models. Chapman and Hall, New York.
Hedley SL, Buckland ST (2004) Spatial models for line transect sampling. Journal of Agricultural, Biological, and Environmental Statistics, 9(2):181-199   
Elenga G, Bonenfant C, Péron G (2020) Distance sampling of duikers in the rainforest: Dealing with transect avoidance. PLOS ONE, 15(10): e0240049   
Miller DL, Burt ML, Rexstad EA, Thomas L (2013) Spatial models for distance sampling data: recent developments and future directions. Methods in Ecology and Evolution, 4(11):1001-1010

