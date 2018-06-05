# DoseExpansion

This project includes five R scripts and one stan file, described as follows in the order that they are called

*decEfficacy.R* This function sources three helper scripts (*genParams.R*, *Functions.R*, and *titesimfunctions_phil_V2.R*) which provide the suite of functions necessary to implement the simulation study in the AOC paper as well as a list called `arglist` as long as the number of unique simulation studies. By default, this is 60 (10 true toxicity and efficacy curves, two choices of skeleton, and three DEC/sample size configuations). Set `array_id` to be any number between 1 and 60, and this script will run a complete set of `niter` replications of this setting, saving the results in an `R` binary file called [array_id]`.R`
