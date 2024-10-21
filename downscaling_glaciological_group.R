# Project: Andean glacier-climatic interactions under solar radiation modification geoengineering [SRMG] //
# Authors: Alfonso Fernandez; Francisco Manquehual-Cheuque; Marcelo Somos-Valenzuela                    //
# Contact: alfernandez@udec.cl; f.manquehual01@ufromail.cl; marcelo.somos@ufrontera.cl                 //
# Last modification: 21-10-2024                                                                       //
# /////////////////////////////////////////////////////////////////////////////////////////////////////


################### MODIFY BY USER ###################

# HOW TO RUN: Rscript downscaling_glaciological_group.R --glaciological_group=group_FA --scenario=historical,ssp245 --target_variable=tas --type_tas=tas --type_prediction=quantity --apply_cv=FALSE --folds=6 --by_season=FALSE --export_netcdf=TRUE --to_celsius=TRUE

# Directories
dir_layers <- '/home/user/layers'
dir_obs <- '/home/user/input/era5'
dir_scenarios <- '/home/user/input/ssp245'
dir_binomial <- '/home/user/output/binomial'
dir_quantity <- '/home/user/output/quantity'

################### END ###################




################### AVOID MAKING CHANGES ###################

# Previous step ----

# Try the following if Java runs out of memory:
options(java.parameters = "-Xmx8000m")

message('------------------------- Loading libraries -------------------------')

suppressPackageStartupMessages(library(loadeR))
suppressPackageStartupMessages(library(loadeR.2nc))
suppressPackageStartupMessages(library(downscaleR))
suppressPackageStartupMessages(library(convertR))
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(rgeos))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

# end ----




# Arguments ----

option_list = list(

  make_option('--glaciological_group', type='character', default='group_FA',
              help="Shapefile layer of the group.", metavar="character"),

  make_option(opt_str='--scenario', type='character', default='all',
              help="Scenario to apply downscaling.
                    Options: historical, ssp245, ssp585, g6solar, all", metavar="character"),

  make_option(opt_str='--target_variable', type='character', default='pr',
              help="Options: pr or tas (if 'tas', remember to define '--temperature_type')", metavar="character"),

  make_option(opt_str='--type_tas', type='character', default='tas',
              help="Options: tas, tasmin, or tasmax", metavar="character"),

  make_option(opt_str='--type_prediction', type='character', default='amount',
              help="Options: quantity or binary", metavar="character"),

  make_option(opt_str='--apply_cv', type='logical', default=FALSE,
              help="Apply cross-validation.
                    Options: TRUE or FALSE", metavar="logical"),

  make_option(opt_str='--folds', type='double', default=6,
              help="Options: Default is 6,
                    following the work of Araya-Osses et al. (2020) (for 36 years)", metavar="double"),

  make_option(opt_str='--by_season', type='logical', default=FALSE,
              help="Options: TRUE or FALSE", metavar="logical"),

  make_option(opt_str='--export_netcdf', type='logical', default=TRUE,
              help="Options: TRUE or FALSE", metavar="logical"),

  make_option(opt_str='--to_celsius', type='logical', default=FALSE,
              help="Kelvin to Celsius.
                    Options: TRUE or FALSE", metavar="logical")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# end ---




# Reading arguments ----

message('------------------------- Reading arguments -------------------------')

glaciological_group <- opt$glaciological_group
scenario <- unlist(strsplit(opt$scenario, split = ","))
target_variable <- opt$target_variable
type_tas <- opt$type_tas
type_prediction <- opt$type_prediction
apply_cv <- opt$apply_cv
folds <- opt$folds
by_season <- opt$by_season
export_netcdf <- opt$export_netcdf
to_celsius <- opt$to_celsius

# Read shp
setwd(dir_layers)
domain <- readOGR('.', glaciological_group, verbose = FALSE)
latitude <- domain@bbox[2,]
longitude <- domain@bbox[1,]


# Setting time to downscaling
if(by_season){
  n_season <- 4
  season <- list(c(1,2,12),
                     c(3,4,5),
                     c(6,7,8),
                     c(9,10,11))
  } else if(!by_season){
    n_season <- 1
    season <- list(1:12)
    } else(stop("Object 'by_season' must be TRUE or FALSE"))


# Part of name file output
name_data_training <- 'ERA5-Historical'

# end ---




# Downscaling process ----

for (j in 1:n_season) {
  months_j <- season[[j]]

  if(j==1 & by_season==FALSE){season_j <- 'annual'
  } else if(j==1 & by_season==TRUE){season_j <- 'DJF'
  } else if(j==2){season_j <- 'MAM'
  } else if(j==3){season_j <- 'JJA'
  } else if(j==4){season_j <- 'SON'}

  message('------------------------- Season: ', season_j,' -------------------------')

  for (i in 1:length(scenario)) {

    # Preparing parameters  ----

    message('------------------------- Type of prediction: ',type_prediction,' -------------------------')
    scenario_i <- scenario[i]

    message('------------------------- Scenario: ',scenario_i,' -------------------------')
    message('------------------------- Checking data -------------------------')

    if(target_variable!='tas' & target_variable!='pr'){stop("Parameter 'target_variable' not accepted")}
    if(type_prediction!='quantity' & type_prediction!='binary'){stop("Parameter 'type_prediction' not accepted")}
    if(scenario_i!='historical' & scenario_i!='ssp245' & scenario_i!='ssp585' & scenario_i!='g6solar'){stop("Parametro 'scenario_i' no aceptado")}

    if(i==1){
      if(target_variable=='pr'){threshold <- 0.1}
      years <- 1979:2014
      months_merged <- str_c(months_j,collapse = '_')
    }

    if(scenario_i=='historical'){years.historical <- 1979:2014
    } else if(scenario_i=='ssp245' | scenario_i=='ssp585'){years.historical <- 2015:2100
    } else if(scenario_i=='g6solar'){years.historical <- 2015:2100 # 2020:2100
    } else(stop("Parameter 'scenario_i' not accepted"))

    period.training.downscaling <- paste0('_', years[1],'_' , years[length(years)])
    period.application.downscaling <- paste0('_', years.historical[1],'_' , years.historical[length(years.historical)])

    period.training <- paste0('training period ', years[1], '-' , years[length(years)])
    if(scenario_i=='g6solar'){period.downscaling <- paste0('and downscaling period ', years.historical[1], '-' , years.historical[length(years.historical)])
    } else if(scenario_i!='g6solar'){period.downscaling <- paste0('and downscaling period ', years[length(1)], '-' , years.historical[length(years.historical)])}

    # end ---




    # Reading observational data ----

    if(i==1){

      message('------------------------- Reading observed data -------------------------')
      setwd(dir_obs)

      # huss 700 hPa
      if(target_variable=='pr'){
        message('Reading HUSS 700 HPA...')
        huss_700 <- 'specific_humidity_700hPa_day.nc'
        huss_700 <- loadGridData(dataset = huss_700,
                                 var = "q",
                                 lonLim = longitude,
                                 latLim= latitude,
                                 season= months_j,
                                 years = years)

        attr(huss_700$Variable, 'units') <- 'kg.kg-1'
        huss_700$Variable$varName <- 'huss_700'
      }


      # tasmin
      if(target_variable=='tas' & type_tas=='tasmin'){
        message('Reading TASMIN ...')
        tasmin <- '2m_minimum_temperature_day.nc'
        tasmin <- loadGridData(dataset = tasmin,
                               var = "t2m",
                               lonLim = longitude,
                               latLim= latitude,
                               season= months_j,
                               years = years)

        attr(tasmin$Variable, 'units') <- 'K'
        tasmin$Variable$varName <- 'tasmin'
        if(to_celsius){tasmin <- udConvertGrid(tasmin, 'celsius')}
      }


      # tasmax
      if(target_variable=='tas' & type_tas=='tasmax'){
        message('Reading TASMAX...')
        tasmax <- '2m_maximum_temperature_day.nc'
        tasmax <- loadGridData(dataset = tasmax,
                               var = "t2m",
                               lonLim = longitude,
                               latLim= latitude,
                               season= months_j,
                               years = years)

        attr(tasmax$Variable, 'units') <- 'K'
        tasmax$Variable$varName <- 'tasmax'
        if(to_celsius){tasmax <- udConvertGrid(tasmax, 'celsius')}
      }


      # tas
      if(target_variable=='pr' | target_variable=='tas'){
        message('Reading TAS ...')
        tas <- '2m_temperature_day.nc'
        tas <- loadGridData(dataset = tas,
                            var = "t2m",
                            lonLim = longitude,
                            latLim= latitude,
                            season= months_j,
                            years = years)

        attr(tas$Variable, 'units') <- 'K'
        tas$Variable$varName <- 'tas'
        if(to_celsius){tas <- udConvertGrid(tas, 'celsius')}
      }


      # t700
      if(target_variable=='tas'){
        message('Reading T700 ...')
        t700 <- 'temperature_700hPa_day.nc'
        t700 <- loadGridData(dataset = t700,
                            var = "t",
                            lonLim = longitude,
                            latLim= latitude,
                            season= months_j,
                            years = years)

        attr(t700$Variable, 'units') <- 'K'
        t700$Variable$varName <- 'ta_700'
        if(to_celsius){t700 <- udConvertGrid(t700, 'celsius')}
      }


      # t850
      if(target_variable=='tas'){
        message('Reading T850 ...')
        t850 <- 'temperature_850hPa_day.nc'
        t850 <- loadGridData(dataset = t850,
                             var = "t",
                             lonLim = longitude,
                             latLim= latitude,
                             season= months_j,
                             years = years)

        attr(t850$Variable, 'units') <- 'K'
        t850$Variable$varName <- 'ta_850'
        if(to_celsius){t850 <- udConvertGrid(t850, 'celsius')}
      }


      # pr
      if(target_variable=='pr'){

        message('Reading PR ...')
        pr <- 'total_precipitation_day.nc'
        pr <- loadGridData(dataset = pr,
                           var = "tp",
                           lonLim = longitude,
                           latLim= latitude,
                           season= months_j,
                           years = years)

        attr(pr$Variable, 'units') <- 'm'
        pr <- udConvertGrid(pr, 'mm')
      }


      # zg 250 hPa
      if(target_variable=='tas'){
        message('Reading ZG 250 HPA ...')
        zg_250 <- 'zg_250hPa_day.nc'
        zg_250 <- loadGridData(dataset = zg_250,
                               var = "z",
                               lonLim = longitude,
                               latLim= latitude-0.25,
                               season= months_j,
                               years = years)

        zg_250$Variable$varName <- 'zg_250hPa'
      }


      # zg 500 hPa
      if(target_variable=='tas' | target_variable=='pr'){
        message('Reading ZG 500 HPA ...')
        zg_500 <- 'zg_500hPa_day.nc'
        zg_500 <- loadGridData(dataset = zg_500,
                               var = "z",
                               lonLim = longitude,
                               latLim= latitude-0.25,
                               season= months_j,
                               years = years)

        zg_500$Variable$varName <- 'zg_500hPa'
      }


      # zg 850 hPa
      if(target_variable=='pr'){
        message('Reading ZG 850 HPA ...')
        zg_850 <- 'zg_850hPa_day.nc'
        zg_850 <- loadGridData(dataset = zg_850,
                               var = "z",
                               lonLim = longitude,
                               latLim= latitude-0.25,
                               season= months_j,
                               years = years)

        zg_850$Variable$varName <- 'zg_850hPa'
      }

      # end ---




      # Observed data preparation ----

      message('------------------------- Preparing data -------------------------')

      if(target_variable=='tas'){t <- eval(parse(text = type_tas))}

      if(target_variable=='pr'){

        predictors <- makeMultiGrid(huss_700, tas, zg_500, zg_850, skip.temporal.check = FALSE)
        y_var <- binaryGrid(pr, condition = "GE", threshold = threshold, partial = TRUE)
        y_bin <- binaryGrid(pr, condition = "GE", threshold = threshold)

      } else if(target_variable=='tas'){

        predictors <- makeMultiGrid(t, t700, t850, zg_250, zg_500, skip.temporal.check = FALSE)

      }

      spatial.pars.M1 <- list(which.combine = getVarNames(predictors), v.exp = .95, rot = FALSE)

      # Before performing PCA, all the variables have to be standardised, since PCA is
      # sensitive to data that has not been centered. That is why one of the preprocessing
      # steps in PCA is to normalise data so that it has μ=0 and σ=1
      # Source: https://rpubs.com/esobolewska/pcr-step-by-step

      # Standardization
      predictors_stand <- scaleGrid(predictors, type = "standardize")
      print(getVarNames(predictors_stand))
    }

    # end ---




    # Reading simulated data ----

    if(!apply_cv){
      message('------------------------- Reading scenario data -------------------------')

      latitude <- domain@bbox[2,]
      longitude <- domain@bbox[1,]

      if(target_variable=='pr'){resolution.reference <- pr}

      if(target_variable=='tas'){
        t <- eval(parse(text = type_tas))
        resolution.reference <- t}


      # huss 700 hPa
      setwd(dir_scenarios)

      if(target_variable=='pr'){

        message('Reading HUSS 700 HPA ...')
        huss_700.historical <- 'huss_700hPa_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        huss_700.historical <- loadGridData(dataset = huss_700.historical,
                                            var = 'hus',
                                            lonLim = longitude+c(-1,1),
                                            latLim= latitude,
                                            season= months_j,
                                            years = years.historical)

        huss_700.historical$Variable$varName <- 'huss_700'
        huss_700.historical <- interpGrid(huss_700.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }


      # tasmin
      if(target_variable=='tas' & type_tas=='tasmin'){

        message('Reading TASMIN ...')
        tasmin.historical <- 'tasmin_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        tasmin.historical <- loadGridData(dataset = tasmin.historical,
                                          var = 'tasmin',
                                          lonLim = longitude+c(-1,1),
                                          latLim= latitude,
                                          season= months_j,
                                          years = years.historical)

        if(to_celsius){tasmin.historical <- udConvertGrid(tasmin.historical, 'celsius')}
        tasmin.historical$Variable$varName <- 'tasmin'
        tasmin.historical <- interpGrid(tasmin.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }


      # tasmax
      if(target_variable=='tas' & type_tas=='tasmax'){

        message('Reading TASMAX ...')
        tasmax.historical <- 'tasmax_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        tasmax.historical <- loadGridData(dataset = tasmax.historical,
                                          var = 'tasmax',
                                          lonLim = longitude+c(-1,1),
                                          latLim= latitude,
                                          season= months_j,
                                          years = years.historical)

        if(to_celsius){tasmax.historical <- udConvertGrid(tasmax.historical, 'celsius')}
        tasmax.historical$Variable$varName <- 'tasmax'
        tasmax.historical <- interpGrid(tasmax.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }


      # tas
      if(target_variable=='pr' | target_variable=='tas'){

        message('Reading TAS ...')
        tas.historical <- 'tas_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        tas.historical <- loadGridData(dataset = tas.historical,
                                       var = 'tas',
                                       lonLim = longitude+c(-1,1),
                                       latLim= latitude,
                                       season= months_j,
                                       years = years.historical)

        if(to_celsius){tas.historical <- udConvertGrid(tas.historical, 'celsius')}
        tas.historical$Variable$varName <- 'tas'
        tas.historical <- interpGrid(tas.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }


      # t700
      if(target_variable=='tas'){

        message('Reading T700 HPA ...')
        t700.historical <- 'ta_700hPa_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        t700.historical <- loadGridData(dataset = t700.historical,
                                       var = 'ta',
                                       lonLim = longitude+c(-1,1),
                                       latLim= latitude+c(-1,1),
                                       season= months_j,
                                       years = years.historical)


        if(to_celsius){t700.historical <- udConvertGrid(t700.historical, 'celsius')}
        t700.historical$Variable$varName <- 'ta_700'
        t700.historical <- interpGrid(t700.historical, new.coordinates = getGrid(resolution.reference), method = "nearest",
                                      force.non.overlapping=TRUE)
      }


      # t850
      if(target_variable=='tas'){

        message('Reading T850 HPA ...')
        t850.historical <- 'ta_850hPa_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        t850.historical <- loadGridData(dataset = t850.historical,
                                        var = 'ta',
                                        lonLim = longitude+c(-1,1),
                                        latLim= latitude,
                                        season= months_j,
                                        years = years.historical)

        if(to_celsius){t850.historical <- udConvertGrid(t850.historical, 'celsius')}
        t850.historical$Variable$varName <- 'ta_850'
        t850.historical <- interpGrid(t850.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }


      # zg 250 hPa
      if(target_variable=='tas'){

        message('Reading ZG 250 HPA ...')
        zg_250.historical <- 'zg_250hPa_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        zg_250.historical <- loadGridData(dataset = zg_250.historical,
                                          var = 'zg',
                                          lonLim = longitude+c(-1,1),
                                          latLim= latitude,
                                          season= months_j,
                                          years = years.historical)

        zg_250.historical$Variable$varName <- 'zg_250hPa'
        zg_250.historical <- interpGrid(zg_250.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }


      # zg 500 hPa
      if(target_variable=='pr' | target_variable=='tas'){

        message('Reading ZG 500 HPA ...')
        zg_500.historical <- 'zg_500hPa_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        zg_500.historical <- loadGridData(dataset = zg_500.historical,
                                          var = 'zg',
                                          lonLim = longitude+c(-1,1),
                                          latLim= latitude,
                                          season= months_j,
                                          years = years.historical)

        zg_500.historical$Variable$varName <- 'zg_500hPa'
        zg_500.historical <- interpGrid(zg_500.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }


      # zg 850 hPa
      if(target_variable=='pr'){

        message('Reading ZG 850 HPA ...')
        zg_850.historical <- 'zg_850hPa_CNRM-ESM2-1_historical_ssp245_1979_2100.nc'
        zg_850.historical <- loadGridData(dataset = zg_850.historical,
                                          var = 'zg',
                                          lonLim = longitude+c(-1,1),
                                          latLim= latitude,
                                          season= months_j,
                                          years = years.historical)

        zg_850.historical$Variable$varName <- 'zg_850hPa'
        zg_850.historical <- interpGrid(zg_850.historical, new.coordinates = getGrid(resolution.reference), method = "nearest")
      }

      # end ---




      # Preparation of simulated data ----

      message('------------------------- Preparing data -------------------------')

      # Collection of predictors
      type_tas.simulated <- paste0(type_tas,'.historical')

      if(target_variable=='pr'){
        predictors.historical0 <- makeMultiGrid(huss_700.historical, tas.historical, zg_500.historical,
                                                 zg_850.historical, skip.temporal.check = FALSE)}
      if(target_variable=='tas'){
        t.historical <- eval(parse(text = type_tas.simulated))
        predictors.historical0 <- makeMultiGrid(t.historical, t700.historical, t850.historical,
                                                 zg_250.historical, zg_500.historical, skip.temporal.check = TRUE)}
      print(getVarNames(predictors.historical0))

      # Scaling values
      predictors.historical <- scaleGrid(predictors.historical0, base = predictors.historical0, ref = predictors,
                                          type = "center", # 'center' to Temperature and 'ratio' for Precipitation
                                          spatial.frame = "gridbox",time.frame = "daily")
      rm(list = 'predictors.historical0')

      predictors.historical <- scaleGrid(predictors.historical, base = predictors,type = "standardize",
                                          skip.season.check = TRUE)
    }

    # end ---




    # Fit-Prediction Binomial ----

    if(type_prediction=='binary'){
      if(target_variable!='pr'){stop('Binary prediction is only made for precipitation')}

      # Fitting
      if(i==1){
        message('------------------------- Fitting model for binary -------------------------')

        if(apply_cv){
          binomial.regression.model <- downscaleCV(x=predictors_stand, y=y_bin, method = "GLM",
                                                    family = binomial(link = "logit"),
                                                    folds = folds, sampling.strategy	= "kfold.random",
                                                    prepareData.args = list(
                                                      "spatial.predictors" = spatial.pars.M1,
                                                      "combined.only" = TRUE))

          historical.binomial <- subsetGrid(binomial.regression.model, var='bin')

        } else if(!apply_cv){
        data.binomial <- prepareData(x=predictors_stand, y=y_bin,
                                     spatial.predictors=spatial.pars.M1, combined.only = TRUE)
        binomial.regression.model <- downscaleTrain(data.binomial, method = "GLM", family = binomial(link = "logit"),
                                                     model.verbose = TRUE)

        } else(message("------------------------- apply_cv' must be TRUE or FALSE -------------------------"))

      } else(message('------------------------- Using fitted model -------------------------'))


      # Prediction
      message('------------------------- Applying downscaling -------------------------')

      if(!apply_cv){
        newdata.binomial <- prepareNewData(predictors.historical, data.binomial)
        if(scenario_i!='historical'){rm(ls=predictors.historical)}

        historical.binomial0 <- downscalePredict(newdata.binomial, binomial.regression.model)
        historical.binomial <- binaryGrid(historical.binomial0, ref.obs = y_bin,
                                          ref.pred = binomial.regression.model$pred)
      }

      # Save outputs
      message('------------------------- Save results -------------------------')

      # Name of output file:
      setwd(dir_binomial)
      if(target_variable=='tas'){output.variable.name <- type_tas} else(output.variable.name <- target_variable)

      if(!by_season){period.application.downscaling2 <- period.application.downscaling
      } else if(by_season==TRUE & months_merged=='12_1_2'){period.application.downscaling2 <- paste0(period.application.downscaling,'_DJF')
      } else if(by_season==TRUE & months_merged=='1_2_12'){period.application.downscaling2 <- paste0(period.application.downscaling,'_DJF')
      } else if(by_season==TRUE & months_merged=='3_4_5'){period.application.downscaling2 <- paste0(period.application.downscaling,'_MAM')
      } else if(by_season==TRUE & months_merged=='6_7_8'){period.application.downscaling2 <- paste0(period.application.downscaling,'_JJA')
      } else if(by_season==TRUE & months_merged=='9_10_11'){period.application.downscaling2 <- paste0(period.application.downscaling,'_SON')
      } else(stop("Parameter 'months_j' not accepted"))

      if(apply_cv){type.downscaling <- '_and_downscaling_binomial_CV'} else(type.downscaling <- '_and_downscaling_binomial')

      fileName <- paste0(output.variable.name,'_training', period.training.downscaling, type.downscaling,
                         period.application.downscaling2, '_with_data_from_', scenario_i, '.nc') ; fileName


      # Including a global attribute:
      globalAttributeList1 <- list("Project" = "Andean glacier-climatic interactions under solar radiation modification geoengineering [SRMG]")
      globalAttributeList2 <- list("Downscaled with" = paste(name_data_training, period.training, period.downscaling))


      # Including variable attributes:
      varAttributeList <- list(var_attr1 = 'mm')


      # Others
      historical.binomial$Variable$varName <- target_variable


      # Create file:
      if(export_netcdf){
        message('Exporting netcdf ...')
        grid2nc(data = historical.binomial,
                NetCDFOutFile = fileName,
                missval = 1e20,
                prec = "float",
                globalAttributes = c(globalAttributeList1, globalAttributeList2),
                varAttributes = varAttributeList)
      } else(message("You decided not to export (export_netcdf=FALSE)"))

    }

    # end ---




    # Fit-Prediction quantity ----

    if(type_prediction=='quantity'){

      if(i==1){
        message('------------------------- Fitting model for quantity -------------------------')

        # Fitting
        if(target_variable=='pr'){

          if(apply_cv){
            gamma.regression.model <- downscaleCV(x=predictors_stand, y=pr, method = "GLM",
                                                   family = Gamma(link = "log"),
                                                   folds = folds, sampling.strategy = "kfold.random",
                                                   condition = "GE", threshold = threshold,
                                                   prepareData.args = list(
                                                     "spatial.predictors" = spatial.pars.M1,
                                                     "combined.only" = TRUE))
          } else if(!apply_cv){
        data.gamma <- prepareData(x=predictors_stand,y=pr,spatial.predictors=spatial.pars.M1, combined.only = TRUE)
        gamma.regression.model <- downscaleTrain(data.gamma, method = "GLM", family = Gamma(link = "log"),
                                                  condition = "GE", threshold = threshold)
          } else(message("------------------------- 'apply_cv' must be TRUE or FALSE -------------------------"))

        } else if(target_variable=='tas'){

          if(type_tas=='tas'){y_obs <- tas
            } else if(type_tas=='tasmin'){y_obs <- tasmin
            } else if(type_tas=='tasmax'){y_obs <- tasmax
            } else(stop("problems with the 'type_tas' object"))

          if(apply_cv){
            gamma.regression.model <- downscaleCV(x=predictors_stand, y=y_obs, method = "GLM",
                                                   family = 'gaussian',
                                                   folds = folds, sampling.strategy	= "kfold.random",
                                                   prepareData.args = list(
                                                     "spatial.predictors" = spatial.pars.M1,
                                                     "combined.only" = TRUE))
          } else if(!apply_cv){
            data.gamma <- prepareData(x=predictors_stand,y=y_obs,spatial.predictors=spatial.pars.M1, combined.only = TRUE)
            gamma.regression.model <- downscaleTrain(data.gamma, method = "GLM", family = 'gaussian')
          } else(message("------------------------- 'apply_cv' must be TRUE or FALSE -------------------------"))


        } else(stop("The 'target_variable' object must be 'pr' or 'tas'"))

      } else(message('------------------------- Using fitted model for quantity -------------------------'))


      # Prediction
      message('------------------------- Applying downscaling -------------------------')

      if(!apply_cv){

        message('Applying downscaling ...')
        newdata.gamma <- prepareNewData(predictors.historical, data.gamma)
        historical.gamma <- downscalePredict(newdata.gamma, gamma.regression.model)
        rm(list='newdata.gamma')

        } else(historical.gamma <- gamma.regression.model)


      # Correcting dates in temperature
      if(target_variable=='tas'){
        historical.gamma$Dates$start <- as.Date(historical.gamma$Dates$start) + 1
        historical.gamma$Dates$end <- as.Date(historical.gamma$Dates$end) + 1
        }


      # Save outputs
      message('------------------------- Save results -------------------------')

      # Name of output
      setwd(dir_quantity)
      if(target_variable=='tas'){output.variable.name <- type_tas} else(output.variable.name <- target_variable)

      if(!by_season){period.application.downscaling2 <- period.application.downscaling
      } else if(by_season==TRUE & months_merged=='12_1_2'){period.application.downscaling2 <- paste0(period.application.downscaling,'_DJF')
      } else if(by_season==TRUE & months_merged=='1_2_12'){period.application.downscaling2 <- paste0(period.application.downscaling,'_DJF')
      } else if(by_season==TRUE & months_merged=='3_4_5'){period.application.downscaling2 <- paste0(period.application.downscaling,'_MAM')
      } else if(by_season==TRUE & months_merged=='6_7_8'){period.application.downscaling2 <- paste0(period.application.downscaling,'_JJA')
      } else if(by_season==TRUE & months_merged=='9_10_11'){period.application.downscaling2 <- paste0(period.application.downscaling,'_SON')
      } else(stop("Parameter 'months_j' not accepted"))

      if(apply_cv){type.downscaling <- '_and_downscaling_quantity_CV'} else(type.downscaling <- '_and_downscaling_quantity')

      fileName <- paste0(output.variable.name,'_training', period.training.downscaling, type.downscaling,
                         period.application.downscaling2, '_with_data_from_', scenario_i, '.nc') ; fileName


      # Including a global attribute:
      globalAttributeList1 <- list("Project" = "Andean glacier-climatic interactions under solar radiation modification geoengineering [SRMG]")
      globalAttributeList2 <- list("Downscaled with" = paste(name_data_training, period.training, period.downscaling))


      # Including variable attributes:
      if(target_variable=='pr'){varAttributeList <- list(var_attr1 = 'mm')
      } else(varAttributeList <- list(var_attr1 = 'celsius'))


      # Others
      historical.gamma$Variable$varName <- target_variable


      # Create file:
      if(export_netcdf){
        message('Exporting netcdf ...')
        grid2nc(data = historical.gamma,
                NetCDFOutFile = fileName,
                missval = 1e20,
                prec = "float",
                globalAttributes = c(globalAttributeList1, globalAttributeList2),
                varAttributes = varAttributeList)
      } else(message("You decided not to export (export_netcdf=FALSE)"))

    }

    # end ---

    message('------------------------- DONE! -------------------------')

  }

}

################### END ###################
