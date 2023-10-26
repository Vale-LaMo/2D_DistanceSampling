# 2D_DistanceSampling
 
In distance sampling surveys, the animals might avoid both the transects in the absence of observers, and the observers themselves. To correct for the effect of the behavioral responses of the animals to either the transects or the observers, we can estimate density and abundance using line transect survey data with both the forward and perpendicular distances to the observers (2D distance sampling - R LT2D package, Borchers and Cox 2017), not just the perpendicular distance. This analysis approach was also applied and recommended by Elenga et al. (2020).

Here, we rely on the functions from LT2D package (https://github.com/david-borchers/LT2D), as partly revised and applied in Elenga et al. (2020) (https://github.com/cbonenfant/duikers-abundance). With respect to the latter, we made additional minor changes to the code.

As an example, we apply the code to estimate the density of two species with different behavioral response to the transects/observers - impala and duikers.

## Density Surface Modelling

*(work in progress)*


## Instructions

In order to replicate the analyses, or to apply them to your datasets, fork this repository and then clone it to your machine. This will re-create the working environment, with the main file *2D_DistanceSampling.Rproj* to be opened in RStudio.
The code for the analyses is stored in the *.Rmd* files *00_LT2D* and *01_dsm*.    
The input data are in the *data* folder.
The *functions* folder includes the scripts with the (customized) R functions.   
The *output* folder is meant to store the results of the analyses (mainly fitted models but also datasets produced in the first stage of the analysis and then used for DSMs).   

The list of packages needed for the analyses is reported at the beginnin of each *.Rmd* file. Please make sure to install them before running the code. The LT2D package should be installed via devtools.

