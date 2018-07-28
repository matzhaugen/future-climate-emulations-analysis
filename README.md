# future-climate-emulations-analysis

This repository contain the code necessary to produce the figures to the paper [CITE].


# Reproducing the temperature projections + figures

To run this code, first clone the repo,
```
git clone https://github.com/matzhaugen/future-climate-emulations-analysis.git
cd future-climate-emulations-analysis
```

and download the raw data from Dropbox:
```
curl https://www.dropbox.com/s/a1dtbml8umsovnp/raw_temperatures.RData?dl=1 -O -J -L
curl https://www.dropbox.com/s/842bdwz9fz3s76r/CCSM4_volcanic_1850-2008_prototype1.nc?dl=1 -O -J -L
```

Make sure that R is installed and has the necessary packages given at the top of `emulationHelper.R`.

Run `maps.R` to produce the projected temperatures. This code should really be run on an 8 core machine and is the default in the file. To modify, change the `mc.cores` variable in line 4 of `maps.R`.

Run `figures.R` to produce the figures. This relies on two files 
1. `everything_1463_1463_310_1010.Rdata`: Contains all the raw model coefficients that defines the quantile maps.
2. `maps_1463_1463_310_1010.Rdata`: Contains the temperature projections.

# Reproducing the quantile maps

This requires a CPU with at least 8 cores (preferably 16) with >6Gb/core running R/3.4.3. The code does not work with R versions >3.5 because `chol` was updated in the base packages and stopped working (don't ask me why).

```
Rscript multiMain_jack.R
```
The output will be saved to `everything_1463_1463_310_1010.Rdata`. For a SLURM submission script see `multiMain_jack.sh`. The job takes about 3 hours on a 16 core machine.

# Deep dive into the code

The three central functions in `emulationHelper.R` are 
1. `quantile.map`: estimates a quantile map based on raw temperatures given either from the LENS or the SFK15 model output. 
2. `map.data`: Uses the model coefficients estimated in `quantile.map` to project any data into the future, either into one year or into a range of years. This function can also handle nested output (see point 3 below).
3. `quantile.map.nested`: Produces quantile maps in a nested fashion by first using the output from `quantile.map` and one simulation or a small set of simulations to fine tune the model coefficients with new data. The idea is to first produce an accurate model with a lot of good data, then re-calibrate using data (e.g. temperatures) from a coarser climate model. Note that there are two types of models here, a climate model from which we get raw temperatures and a statistical model which we also call quantile maps. 
