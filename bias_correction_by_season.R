# Project: Andean glacier-climatic interactions under solar radiation modification geoengineering [SRMG] //
# Authors: Alfonso Fernandez; Francisco Manquehual-Cheuque; Marcelo Somos-Valenzuela                    //
# Contact: alfernandez@udec.cl; f.manquehual01@ufromail.cl; marcelo.somos@ufrontera.cl                 //
# Last modification: 21-10-2024                                                                       //
# /////////////////////////////////////////////////////////////////////////////////////////////////////


################### MODIFY BY USER ###################

# HOW TO RUN: Rscript bias_correction_by_season.R --glaciological_group=group_FA --data_train=tas_training_1979_2014_and_downscaling_quantity_1979_2014_with_data_from_historical.nc --data_to_correct=tas_training_1979_2014_and_downscaling_quantity_2015_2100_with_data_from_ssp245.nc --method=isimip --scenario=ssp245 --target_variable=tas --export_csv=TRUE --export_netcdf=TRUE

# Directories
dir_functions <- '/home/user/routines'
dir_layers <- '/home/user/layers'
dir_obs <- '/home/user/input/era5'
dir_scenarios <- '/home/user/output/quantity'
dir_corrected <- '/home/user/output/corrected'

################### END ###################




################### AVOID MAKING CHANGES ###################

# Previous step ----

# Try the following if Java runs out of memory:
options(java.parameters = "-Xmx8000m")

message('------------------------- Loading libraries -------------------------')

suppressPackageStartupMessages(library(loadeR))
suppressPackageStartupMessages(library(loadeR.2nc))
suppressPackageStartupMessages(library(visualizeR))
suppressPackageStartupMessages(library(downscaleR))
suppressPackageStartupMessages(library(convertR))
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(rgeos))
suppressPackageStartupMessages(library(optparse))

# end ----




# Arguments ----

option_list = list(

  make_option('--glaciological_group', type='character', default='group_FA',
              help="Shapefile layer of the group.", metavar="character"),

  make_option(opt_str='--data_train', type='character', default='G5_T_hist.nc',
              help="Data to train the bias correction method.", metavar="character"),

  make_option(opt_str='--data_to_correct', type='character', default='G5_T_g6solar.nc',
              help="New data to apply the bias correction trained.", metavar="character"),

  make_option(opt_str='--method', type='character', default='isimip',
              help="Options: eqm, pqm, gpqm, loci, ptr, qdm, isimip, mva, variance, all.", metavar="character"),

  make_option(opt_str='--scenario', type='character', default='historical',
              help="Options: historical, ssp245, ssp585, g6solar", metavar="character"),

  make_option(opt_str='--target_variable', type='character', default='tas',
              help="Options: pr, tas, tasmin or tasmax", metavar="character"),

  make_option(opt_str='--export_csv', type='logical', default=TRUE,
              help="Options: TRUE or FALSE", metavar="logical"),

  make_option(opt_str='--export_netcdf', type='logical', default=TRUE,
              help="Options: TRUE or FALSE", metavar="logical")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# end ---




# Reading arguments ----

message('------------------------- Reading arguments -------------------------')

# Calling external functions
setwd(dir_functions)
source('function_grid_to_mean_data.R')
source('function_bias_correction_by_season.R')


# Read arguments
setwd(dir_layers)
glaciological_group <- opt$glaciological_group
shp_boundary <- readOGR('.', glaciological_group, verbose = FALSE)
latitude <- shp_boundary@bbox[2,]
longitude <- shp_boundary@bbox[1,]

target.variable <- opt$target_variable
method <- opt$method
scenario <- opt$scenario
data.over <- scenario
export.csv <- opt$export_csv
export.netcdf <- opt$export_netcdf
data_train <- opt$data_train
data_to_correct <- opt$data_to_correct


# Setting more variables
correction.data <- 'ERA5'
window <- NULL
training.years <- 1979:2014

if(data.over=='historical'){years.to.correct <- 1979:2014
} else if(data.over=='ssp245' | data.over=='ssp585'){years.to.correct <- 2015:2100
} else if(data.over=='g6solar'){years.to.correct <- 2020:2100
} else(stop("Parameter 'data.over.i' not accepted"))

if(target.variable=='pr'){output.name <- 'total_precipitation_'
} else if(target.variable=='tas'){output.name <- 'temperature_'
} else if(target.variable=='tasmin'){output.name <- 'minimum_temperature_'
} else if(target.variable=='tasmax'){output.name <- 'maximum_temperature_'
} else(stop("The object 'target.variable' must be 'tas', 'tasmin', or 'tasmax'"))

if(data.over=='g6solar'){downscaling_years <- '2020_2100'
} else if(data.over=='historical'){downscaling_years <- '1979_2014'
} else if(data.over=='ssp245' | data.over=='ssp585'){downscaling_years <- '2015_2100'}

if(target.variable=='pr'){threshold <- 0.1
                            window <- NULL}

training.period <- paste0('training period ', training.years[1], '-' , training.years[length(training.years)])
correction.period <- paste0('and corrected period ', years.to.correct[1], '-' , years.to.correct[length(years.to.correct)])

message('------------------------- Glaciological group: ',glaciological_group,' -------------------------')
message('------------------------- Scenario: ',data.over,'-------------------------')

if(target.variable=='pr'){target.variable.message <- 'Precipitation'
} else if(target.variable=='tas'){target.variable.message <- 'Mean Temperature'
} else if(target.variable=='tasmin'){target.variable.message <- 'Minimum Temperature'
} else if(target.variable=='tasmax'){target.variable.message <- 'Maximum Temperature'}

message('------------------------- Target variable: ',target.variable.message,' -------------------------')

# end ---




# Reading grids ----

# Observed
message('------------------------- Reading observed data -------------------------')

setwd(dir_obs)
if(target.variable=='pr'){nc.name <- 'total_precipitation_day.nc'
} else if(target.variable=='tas'){nc.name <- '2m_temperature_day.nc'
} else if(target.variable=='tasmin'){nc.name <- '2m_minimum_temperature_day.nc'
} else if(target.variable=='tasmax'){nc.name <- '2m_maximum_temperature_day.nc'
} else(stop("Variable provided to 'target.variable' not accepted"))

if(target.variable=='pr'){var.i <- 'tp'
} else if(target.variable=='tas' | target.variable=='tasmin' | target.variable=='tasmax'){var.i <- 't2m'}

# Because the observed data starts on 01-01-1979 and the simulated downscaling on 02-01-1979
if(method=='isimip'){training.years <- 1980:2014}

observed <- loadGridData(dataset = nc.name,
                         var = var.i,
                         years = training.years,
                         lonLim = longitude,
                         latLim = latitude)

if(target.variable=='pr'){
  attr(observed$Variable, 'units') <- 'm'
  observed <- udConvertGrid(observed, 'mm')
} else if(target.variable=='tas' | target.variable=='tasmin' | target.variable=='tasmax'){
  attr(observed$Variable, 'units') <- 'K'
  observed <- udConvertGrid(observed, 'celsius')}


# Simulated
message('------------------------- Reading simulated data for training -------------------------')

setwd(dir_scenarios)

if(target.variable=='pr'){
  var.i <- 'pr'
} else if(target.variable=='tas' | target.variable=='tasmin' | target.variable=='tasmax'){
  var.i <- 'tas'}

simulated <- loadGridData(dataset = data_train, var = var.i, years = training.years)

if(target.variable=='pr' & data.over=='historical' & glaciological_group=='group_IT' | glaciological_group=='group_OT1' | glaciological_group=='group_OT2'){
  message('Correcting dates ...')
  simulated$Dates$start <- paste(as.Date(simulated$Dates$start)-1, '00:00:00 GMT')
  simulated$Dates$end <- paste(as.Date(simulated$Dates$end)-1, '00:00:00 GMT')
}


# Simulated newdata
message('------------------------- Reading simulated data to be corrected -------------------------')

simulated.newdata <- loadGridData(dataset = data_to_correct, var = var.i, years = years.to.correct)

if(target.variable=='pr' & data.over=='historical' & glaciological_group=='group_IT' | glaciological_group=='group_OT1' | glaciological_group=='group_OT2'){
  message('Correcting dates ...')
  simulated.newdata$Dates$start <- paste(as.Date(simulated.newdata$Dates$start)-1, '00:00:00 GMT')
  simulated.newdata$Dates$end <- paste(as.Date(simulated.newdata$Dates$end)-1, '00:00:00 GMT')
}

# end ---




# Bias correction ----

message('------------------------- Starting bias correction -------------------------')

# Must be for the 4 seasons! Source: Manzanas et al 2017
# DJF, MAM, JJA, SON (Araya-Osses et al., 2020)

if(target.variable=='pr'){

  if(method=='eqm' | method=='all'){
    message('------------------------- Method EQM -------------------------')
    eqm <- bias_correction_by_season(method = 'eqm', threshold = threshold, window = window, precipitation = TRUE,
                                     x=simulated, y=observed, y.newdata = simulated.newdata,
                                     training.years = training.years, correction.years = years.to.correct)
  }

  if(method=='pqm' | method=='all'){
    message('------------------------- Method PQM -------------------------')
    pqm <- bias_correction_by_season(method = 'pqm', threshold = threshold, window = window, precipitation = TRUE,
                                     x=simulated, y=observed, y.newdata = simulated.newdata,
                                     training.years = training.years, correction.years = years.to.correct)
  }

  if(method=='gpqm' | method=='all'){
    message('------------------------- Method GPQM -------------------------')
    gpqm <- bias_correction_by_season(method = 'gpqm', threshold = threshold, window = window, precipitation = TRUE,
                                      x=simulated, y=observed, y.newdata = simulated.newdata,
                                      training.years = training.years, correction.years = years.to.correct)
  }

  if(method=='loci' | method=='all'){
    message('------------------------- Method LOCI -------------------------')
    loci <- bias_correction_by_season(method = 'loci', threshold = threshold, window = window, precipitation = TRUE,
                                      x=simulated, y=observed, y.newdata = simulated.newdata,
                                      training.years = training.years, correction.years = years.to.correct)
  }

  if(method=='ptr' | method=='all'){
    message('------------------------- Method PTR -------------------------')
    ptr <- bias_correction_by_season(method = 'ptr', threshold = threshold, window = window, precipitation = TRUE,
                                     x=simulated, y=observed, y.newdata = simulated.newdata,
                                     training.years = training.years, correction.years = years.to.correct)
  }

  if(method=='qdm' | method=='all'){
    message('------------------------- Method QDM -------------------------')
    qdm <- bias_correction_by_season(method = 'qdm', threshold = threshold, window = window, precipitation = TRUE,
                                     x=simulated, y=observed, y.newdata = simulated.newdata,
                                     training.years = training.years, correction.years = years.to.correct)
  }

  if(method=='isimip' | method=='all'){
    message('------------------------- Method ISIMIP -------------------------')
    isimip <- bias_correction_by_season(method = 'isimip', threshold = threshold, window = window, precipitation = TRUE,
                                        x=simulated, y=observed, y.newdata = simulated.newdata,
                                        training.years = training.years, correction.years = years.to.correct)
  }

} else if(target.variable=='tas' | target.variable=='tasmin' | target.variable=='tasmax'){

  # Methods DELTA, SCALING and DQM have issues when applied by season and annually (tasmin and tasmax)

  if(method=='eqm' | method=='all'){
    message('------------------------- Method EQM -------------------------')
    eqm <- bias_correction_by_season(method = 'eqm', window = window, precipitation=FALSE,
                                     x=simulated, y=observed, y.newdata = simulated.newdata,
                                     training.years = training.years, correction.years = years.to.correct)}

  if(method=='qdm' | method=='all'){
    message('------------------------- Method QDM -------------------------')
    qdm <- bias_correction_by_season(method = 'qdm', window = window, precipitation=FALSE,
                                     x=simulated, y=observed, y.newdata = simulated.newdata,
                                     training.years = training.years, correction.years = years.to.correct)}

  if(method=='isimip' | method=='all'){
    message('------------------------- Method ISIMIP -------------------------')
    isimip <- bias_correction_by_season(method = 'isimip', window = window, precipitation=FALSE, isimip.type = 'additive',
                                        x=simulated, y=observed, y.newdata = simulated.newdata,
                                        training.years = training.years, correction.years = years.to.correct)}

  if(method=='mva' | method=='all'){
    message('------------------------- Method MVA -------------------------')
    mva <- bias_correction_by_season(method = 'mva', window = window, precipitation=FALSE,
                                     x=simulated, y=observed, y.newdata = simulated.newdata,
                                     training.years = training.years, correction.years = years.to.correct)}

  if(method=='variance' | method=='all'){
    message('------------------------- Method VARIANCE -------------------------')
    variance <- bias_correction_by_season(method = 'variance', window = window, precipitation=FALSE,
                                          x=simulated, y=observed, y.newdata = simulated.newdata,
                                          training.years = training.years, correction.years = years.to.correct)}
}

# end ---




# Grids to data frames ----

if(export.csv){
  message('------------------------- Corrected grids to data frames -------------------------')

  if(method=='eqm' | method=='all'){
    df.eqm <- grid_to_mean_data(eqm)}
  if(method=='pqm' | method=='all'){
    df.pqm <- grid_to_mean_data(pqm)}
  if(method=='gpqm' | method=='all'){
    df.gpqm <- grid_to_mean_data(gpqm)}
  if(method=='loci' | method=='all'){
    df.loci <- grid_to_mean_data(loci)}
  if(method=='ptr' | method=='all'){
    df.ptr <- grid_to_mean_data(ptr)}
  if(method=='qdm' | method=='all'){
    df.qdm <- grid_to_mean_data(qdm)}
  if(method=='isimip' | method=='all'){
    df.isimip <- grid_to_mean_data(isimip)}
  if(method=='mva' | method=='all'){
    df.mva <- grid_to_mean_data(mva)}
  if(method=='variance' | method=='all'){
    df.variance <- grid_to_mean_data(variance)}


  # Saving files
  message('------------------------- Saving files (csv) -------------------------')

  setwd(dir_corrected)
  output_file_name_part_2 <- paste0(output.name, data.over, '_corrected')
  output_file_name <- paste(glaciological_group, output_file_name_part_2, sep = '_') ; output_file_name

  if(method=='eqm' | method=='all'){
    write.csv(df.eqm, paste0(output_file_name, '_seasonal_EQM.csv'), row.names = FALSE)}
  if(method=='pqm' | method=='all'){
    write.csv(df.pqm, paste0(output_file_name, '_seasonal_PQM.csv'), row.names = FALSE)}
  if(method=='gpqm' | method=='all'){
    write.csv(df.gpqm, paste0(output_file_name, '_seasonal_GPQM.csv'), row.names = FALSE)}
  if(method=='loci' | method=='all'){
    write.csv(df.loci, paste0(output_file_name, '_seasonal_LOCI.csv'), row.names = FALSE)}
  if(method=='ptr' | method=='all'){
    write.csv(df.ptr, paste0(output_file_name, '_seasonal_PTR.csv'), row.names = FALSE)}
  if(method=='qdm' | method=='all'){
    write.csv(df.qdm, paste0(output_file_name, '_seasonal_QDM.csv'), row.names = FALSE)}
  if(method=='isimip' | method=='all'){
    write.csv(df.isimip, paste0(output_file_name, '_seasonal_ISIMIP.csv'), row.names = FALSE)}
  if(method=='mva' | method=='all'){
    write.csv(df.mva, paste0(output_file_name, '_seasonal_MVA.csv'), row.names = FALSE)}
  if(method=='variance' | method=='all'){
    write.csv(df.variance, paste0(output_file_name, '_seasonal_VARIANCE.csv'), row.names = FALSE)}
}

# end ---




# Saving corrected grids ----

if(export.netcdf){

  message('Exporting netcdf ...')

  setwd(dir_corrected)
  output_file_name_part_2 <- paste0(target.variable,'_', data.over, '_corrected')
  output_file_name <- paste(glaciological_group, output_file_name_part_2, sep = '_') ; output_file_name


  # Name of output file:
  if(method=='eqm' | method=='all'){
    fileName_EQM <- paste0(output_file_name, '_seasonal_EQM.nc4')}
  if(method=='pqm' | method=='all'){
    fileName_PQM <- paste0(output_file_name, '_seasonal_PQM.nc4')}
  if(method=='gpqm' | method=='all'){
    fileName_GPQM <- paste0(output_file_name, '_seasonal_GPQM.nc4')}
  if(method=='loci' | method=='all'){
    fileName_LOCI <- paste0(output_file_name, '_seasonal_LOCI.nc4')}
  if(method=='ptr' | method=='all'){
    fileName_PTR <- paste0(output_file_name, '_seasonal_PTR.nc4')}
  if(method=='qdm' | method=='all'){
    fileName_QDM <- paste0(output_file_name, '_seasonal_QDM.nc4')}
  if(method=='isimip' | method=='all'){
    fileName_ISIMIP <- paste0(output_file_name, '_seasonal_ISIMIP.nc4')}
  if(method=='mva' | method=='all'){
    fileName_MVA <- paste0(output_file_name, '_seasonal_MVA.nc4')}
  if(method=='variance' | method=='all'){
    fileName_VARIANCE <- paste0(output_file_name, '_seasonal_VARIANCE.nc4')}


  # Including a global attribute:
  globalAttributeList1 <- list("Project" = "Andean glacier-climatic interactions under solar radiation modification geoengineering [SRMG]")
  globalAttributeList2 <- list("Corrected with" = paste(correction.data, training.period, correction.period))


  # Including variable attributes:
  if(target.variable=='pr'){
    varAttributeList <- list(var_attr1 = "Precipitation (mm)")
  } else if(target.variable=='tas' | target.variable=='tasmin' | target.variable=='tasmax'){
    varAttributeList <- list(var_attr1 = "Temperature (celsius)")
  }


  # Others
  description <- 'downscaled and corrected'

  if(target.variable=='pr'){

    if(method=='eqm' | method=='all'){
      attr(eqm$Variable, 'units') <- 'mm'
      attr(eqm, 'description') <- description}
    if(method=='pqm' | method=='all'){
      attr(pqm$Variable, 'units') <- 'mm'
      attr(pqm, 'description') <- description}
    if(method=='gpqm' | method=='all'){
      attr(gpqm$Variable, 'units') <- 'mm'
      attr(gpqm, 'description') <- description}
    if(method=='loci' | method=='all'){
      attr(loci$Variable, 'units') <- 'mm'
      attr(loci, 'description') <- description}
    if(method=='ptr' | method=='all'){
      attr(ptr$Variable, 'units') <- 'mm'
      attr(ptr, 'description') <- description}
    if(method=='qdm' | method=='all'){
      attr(qdm$Variable, 'units') <- 'mm'
      attr(qdm, 'description') <- description}
    if(method=='isimip' | method=='all'){
      attr(isimip$Variable, 'units') <- 'mm'
      attr(isimip, 'description') <- description}

  } else if(target.variable=='tas' | target.variable=='tasmin' | target.variable=='tasmax'){

    if(method=='eqm' | method=='all'){
      attr(eqm$Variable, 'units') <- 'celsius'
      attr(eqm, 'description') <- description}
    if(method=='qdm' | method=='all'){
      attr(qdm$Variable, 'units') <- 'celsius'
      attr(qdm, 'description') <- description}
    if(method=='isimip' | method=='all'){
      attr(isimip$Variable, 'units') <- 'celsius'
      attr(isimip, 'description') <- description}
    if(method=='mva' | method=='all'){
      attr(mva$Variable, 'units') <- 'celsius'
      attr(mva, 'description') <- description}
    if(method=='variance' | method=='all'){
      attr(variance$Variable, 'units') <- 'celsius'
      attr(variance, 'description') <- description}
  }


  # Create file:
  if(method=='eqm' | method=='all'){
    grid2nc(data = eqm,
            NetCDFOutFile = fileName_EQM,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='pqm' | method=='all'){
    grid2nc(data = pqm,
            NetCDFOutFile = fileName_PQM,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='gpqm' | method=='all'){
    grid2nc(data = gpqm,
            NetCDFOutFile = fileName_GPQM,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='loci' | method=='all'){
    grid2nc(data = loci,
            NetCDFOutFile = fileName_LOCI,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='ptr' | method=='all'){
    grid2nc(data = ptr,
            NetCDFOutFile = fileName_PTR,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='qdm' | method=='all'){
    grid2nc(data = qdm,
            NetCDFOutFile = fileName_QDM,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='isimip' | method=='all'){
    grid2nc(data = isimip,
            NetCDFOutFile = fileName_ISIMIP,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='mva' | method=='all'){
    grid2nc(data = mva,
            NetCDFOutFile = fileName_MVA,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

  if(method=='variance' | method=='all'){
    grid2nc(data = variance,
            NetCDFOutFile = fileName_VARIANCE,
            missval = 1e20,
            prec = "float",
            globalAttributes = c(globalAttributeList1,globalAttributeList2),
            varAttributes = varAttributeList)}

} else(message("You decided not to export (export.netcdf=FALSE)"))

message('------------------------- Done! -------------------------')

################### END ###################
