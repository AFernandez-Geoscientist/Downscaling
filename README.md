# Downscaling and bias-correction scripts 

## Overview

This repository contains R scripts used to downscale and bias-correct climate model data used in the following paper: Fernández, Manquehual-Cheuque & Somos_Valenzuela (2024): Impact of Solar Radiation Management on Andean glacier-wide surface mass balance, npj Climate and Atmospheric Science, doi: 10.1038/s41612-024-00807-x

# Requirements

Before running the statistical downscaling and bias correction routines, the user must install the following R packages:

1. `climate4R`
2. `raster`
3. `rgdal`
4. `rgeos`
5. `stringr`
6. `optparse`

The routines were tested in Ubuntu 18, using the terminal.

---

# Downscaling scripts

To apply statistical downscaling, you must first update the directories (or PATH) that are at the beginning of the routine `downscaling_glaciological_group.R`. Then, ensure you have all the inputs in the corresponding directories.

An example of how to run the routine:

```bash
Rscript downscaling_glaciological_group.R --glaciological_group=group_FA --scenario=historical,ssp245 --target_variable=tas --type_tas=tas --type_prediction=quantity --apply_cv=FALSE --folds=6 --by_season=FALSE --export_netcdf=TRUE --to_celsius=TRUE
```

To see the different options, run the following:

```bash
Rscript downscaling_glaciological_group.R --help
```

# Bias-correction scripts

To apply bias correction, you must first update the directories (or PATH) that are at the beginning of the `bias_correction_by_season.R` routine. Then, ensure you have all the inputs in the corresponding directories.

An example of how to run the routine:

```bash
Rscript bias_correction_by_season.R --glaciological_group=group_FA --data_train=tas_training_data_from_historical.nc --data_to_correct=tas_training_data_from_ssp245.nc --method=isimip --scenario=ssp245 --target_variable=tas --export_csv=TRUE --export_netcdf=TRUE
```

To see the different options, run the following:

```bash
Rscript bias_correction_by_season.R --help
```
