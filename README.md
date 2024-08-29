# dynamic_niches

This code accompanies the paper "Migratory birds modulate niche tradeoffs in rhythm with seasons and life history" by Yanco et al. PNAS http://doi.org/10.1073/pnas.2316827121

The main workflow for this analysis can be run by following the steps oultined in the file `wf_dynamic_niches.sh`.  This produces all anlayses and then primary results plots.  The remaining ancilary plots can be found in the `./src/plots` directory

Environmental annotations rely on the Mosey system - code for that can be found in `src/mosey` and more details can be found at https://benscarlson.github.io/mosey An associated (but nut strictly necessary) conda environment is included in `src/conda` to assist with packages for annotation.  Necessary packages can also be installed independently (and may be an easier choice, especially if you don't already have conda installed and initialized).

`src/funs` and `src/init` contian convenience code called by other scripts.

## Questions?

If you have any questions about this code base or the paper itself, please email:  
Scott Yanco. 
syanco@umich.edu