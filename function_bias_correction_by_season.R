bias_correction_by_season <- function(
           method = 'eqm', threshold = 0.1, x = NULL, y = NULL, y.newdata = NULL,
           window = NULL, training.years = NULL, correction.years = NULL,
           precipitation=FALSE, isimip.type='multiplicative', scaling.type='additive'
    ){

    message(stringr::str_c("Declared 'precipitation=",precipitation,"'"))
    if(method == 'loci' & threshold < 1){threshold <- 1
                                          message('Threshold changed to 1 mm (see LOCI method details)')}




    # DJF ----

    message('---------- Preparing grids for DJF season ----------')

    x.DJF <- subsetGrid(x, season = c(1,2,12))
    y.DJF <- subsetGrid(y, season = c(1,2,12))
    y.newdata.DJF <- subsetGrid(y.newdata, season = c(1,2,12))


    # Correction
    message('---------- Correcting DJF season ----------')

    if(method == 'pqm'){
    DJF.corrected <- biasCorrection(x = x.DJF,
                                    y = y.DJF,
                                    newdata = y.newdata.DJF,
                                    window = window,
                                    fitdistr.args = list(densfun = "gamma"),
                                    precipitation = precipitation,
                                    method = method,
                                    wet.threshold = threshold)

    } else if(method == 'isimip'){
      if(precipitation == FALSE){isimip.type <- 'additive'}
      message(stringr::str_c('Correction type: ',isimip.type))
      DJF.corrected <- isimip(x = x.DJF,
                              y = y.DJF,
                              newdata = y.newdata.DJF,
                              type = isimip.type,
                              threshold = threshold)

    } else if(method == 'scaling' | method == 'dqm'){
      message(stringr::str_c('Correction type: ',scaling.type))
      DJF.corrected <- biasCorrection(x = x.DJF,
                                      y = y.DJF,
                                      newdata = y.newdata.DJF,
                                      method = method,
                                      scaling.type = scaling.type)

    } else(
      DJF.corrected <- biasCorrection(x = x.DJF,
                                      y = y.DJF,
                                      newdata = y.newdata.DJF,
                                      window = window,
                                      precipitation = precipitation,
                                      method = method,
                                      wet.threshold = threshold)
    )

    DJF.corrected.no.NA <- suppressMessages(filterNA(DJF.corrected))

    message('---------- Done ----------')




    # MAM ----

    message('---------- Preparing grids for MAM season ----------')

    # March
    month <- 3

    x.march <- subsetGrid(x, years = training.years, season = 12)
    x.march0 <- subsetGrid(x, years = training.years, season = month)
    x.march$Data <- x.march0$Data

    y.march <- subsetGrid(y, years = training.years, season = 12)
    y.march0 <- subsetGrid(y, years = training.years, season = month)
    y.march$Data <- y.march0$Data

    y.newdata.march <- subsetGrid(y.newdata, years = correction.years, season = 12)
    y.newdata.march0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.march$Data <- y.newdata.march0$Data


    # April
    month <- 4

    x.april <- subsetGrid(x, years = training.years, season = 11)
    x.april0 <- subsetGrid(x, years = training.years, season = month)
    x.april$Data <- x.april0$Data

    y.april <- subsetGrid(y, years = training.years, season = 11)
    y.april0 <- subsetGrid(y, years = training.years, season = month)
    y.april$Data <- y.april0$Data

    y.newdata.april <- subsetGrid(y.newdata, years = correction.years, season = 11)
    y.newdata.april0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.april$Data <- y.newdata.april0$Data


    # May
    month <- 5

    x.may <- subsetGrid(x, years = training.years, season = 1)
    x.may0 <- subsetGrid(x, years = training.years, season = month)
    x.may$Data <- x.may0$Data

    y.may <- subsetGrid(y, years = training.years, season = 1)
    y.may0 <- subsetGrid(y, years = training.years, season = month)
    y.may$Data <- y.may0$Data

    y.newdata.may <- subsetGrid(y.newdata, years = correction.years, season = 1)
    y.newdata.may0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.may$Data <- y.newdata.may0$Data


    # Aggregate
    x.aggregate <- bindGrid(x.march, x.april, x.may, dimension = 'time')
    y.aggregate <- bindGrid(y.march, y.april, y.may, dimension = 'time')
    y.newdata.aggregate <- bindGrid(y.newdata.march, y.newdata.april, y.newdata.may, dimension = 'time')

    # Correction
    message('---------- Correcting MAM season ----------')

    if(method == 'pqm'){
      MAM.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       window = window,
                                       fitdistr.args = list(densfun = "gamma"),
                                       precipitation = precipitation,
                                       method = method,
                                       wet.threshold = threshold)

    } else if(method == 'isimip'){
      message(stringr::str_c('Correction type: ',isimip.type))
      MAM.corrected0 <- isimip(x = x.aggregate,
                               y = y.aggregate,
                               newdata = y.newdata.aggregate,
                               type = isimip.type,
                               threshold = threshold)

    } else if(method == 'scaling' | method == 'dqm'){
      message(stringr::str_c('Correction type: ',scaling.type))
      MAM.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       method = method,
                                       scaling.type = scaling.type)

    } else(
      MAM.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       window = window,
                                       precipitation = precipitation,
                                       method = method,
                                       wet.threshold = threshold)
    )

    # Union grids
    march.corrected <- subsetGrid(MAM.corrected0, season = 12)
    april.corrected <- subsetGrid(MAM.corrected0, season = 11)
    may.corrected <- subsetGrid(MAM.corrected0, season = 1)

    MAM.corrected <- bindGrid(march.corrected, april.corrected, may.corrected, dimension = 'time')
    MAM.corrected$Dates$start <- gsub('-12-', '-03-', MAM.corrected$Dates$start)
    MAM.corrected$Dates$end <- gsub('-12-', '-03-', MAM.corrected$Dates$end)
    MAM.corrected$Dates$start <- gsub('-11-', '-04-', MAM.corrected$Dates$start)
    MAM.corrected$Dates$end <- gsub('-11-', '-04-', MAM.corrected$Dates$end)
    MAM.corrected$Dates$start <- gsub('-01-', '-05-', MAM.corrected$Dates$start)
    MAM.corrected$Dates$end <- gsub('-01-', '-05-', MAM.corrected$Dates$end)

    MAM.corrected.no.NA <- suppressMessages(filterNA(MAM.corrected))

    message('---------- Done ----------')




    # JJA ----

    message('---------- Preparing grids for JJA season ----------')

    # June
    month <- 6

    x.june <- subsetGrid(x, years = training.years, season = 11)
    x.june0 <- subsetGrid(x, years = training.years, season = month)
    x.june$Data <- x.june0$Data

    y.june <- subsetGrid(y, years = training.years, season = 11)
    y.june0 <- subsetGrid(y, years = training.years, season = month)
    y.june$Data <- y.june0$Data

    y.newdata.june <- subsetGrid(y.newdata, years = correction.years, season = 11)
    y.newdata.june0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.june$Data <- y.newdata.june0$Data


    # July
    month <- 7

    x.july <- subsetGrid(x, years = training.years, season = 12)
    x.july0 <- subsetGrid(x, years = training.years, season = month)
    x.july$Data <- x.july0$Data

    y.july <- subsetGrid(y, years = training.years, season = 12)
    y.july0 <- subsetGrid(y, years = training.years, season = month)
    y.july$Data <- y.july0$Data

    y.newdata.july <- subsetGrid(y.newdata, years = correction.years, season = 12)
    y.newdata.july0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.july$Data <- y.newdata.july0$Data


    # August
    month <- 8

    x.august <- subsetGrid(x, years = training.years, season = 1)
    x.august0 <- subsetGrid(x, years = training.years, season = month)
    x.august$Data <- x.august0$Data

    y.august <- subsetGrid(y, years = training.years, season = 1)
    y.august0 <- subsetGrid(y, years = training.years, season = month)
    y.august$Data <- y.august0$Data

    y.newdata.august <- subsetGrid(y.newdata, years = correction.years, season = 1)
    y.newdata.august0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.august$Data <- y.newdata.august0$Data


    # Aggregate
    x.aggregate <- bindGrid(x.june, x.july, x.august, dimension = 'time')
    y.aggregate <- bindGrid(y.june, y.july, y.august, dimension = 'time')
    y.newdata.aggregate <- bindGrid(y.newdata.june, y.newdata.july, y.newdata.august, dimension = 'time')

    # Correction
    message('---------- Correcting JJA season ----------')

    if(method == 'pqm'){
      JJA.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       window = window,
                                       fitdistr.args = list(densfun = "gamma"),
                                       precipitation = precipitation,
                                       method = method,
                                       wet.threshold = threshold)

    } else if(method == 'isimip'){
      message(stringr::str_c('Correction type: ', isimip.type))
      JJA.corrected0 <- isimip(x = x.aggregate,
                               y = y.aggregate,
                               newdata = y.newdata.aggregate,
                               type = isimip.type,
                               threshold = threshold)

    } else if(method == 'scaling' | method == 'dqm'){
      message(stringr::str_c('Correction type: ', scaling.type))
      JJA.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       method = method,
                                       scaling.type = scaling.type)

    } else(
      JJA.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       window = window,
                                       precipitation = precipitation,
                                       method = method,
                                       wet.threshold = threshold)
    )


    # Union grid
    june.corrected <- subsetGrid(JJA.corrected0, season = 11)

    july.corrected <- subsetGrid(JJA.corrected0, season = 12)

    august.corrected <- subsetGrid(JJA.corrected0, season = 1)

    JJA.corrected <- bindGrid(june.corrected, july.corrected, august.corrected, dimension = 'time')

    JJA.corrected$Dates$start <- gsub('-11-', '-06-', JJA.corrected$Dates$start)
    JJA.corrected$Dates$end <- gsub('-11-', '-06-', JJA.corrected$Dates$end)
    JJA.corrected$Dates$start <- gsub('-12-', '-07-', JJA.corrected$Dates$start)
    JJA.corrected$Dates$end <- gsub('-12-', '-07-', JJA.corrected$Dates$end)
    JJA.corrected$Dates$start <- gsub('-01-', '-08-', JJA.corrected$Dates$start)
    JJA.corrected$Dates$end <- gsub('-01-', '-08-', JJA.corrected$Dates$end)

    JJA.corrected.no.NA <- suppressMessages(filterNA(JJA.corrected))

    message('---------- Done ----------')




    # SON ----

    message('---------- Preparing grids for SON season ----------')

    # September
    month <- 9

    x.september0 <- subsetGrid(x, years = training.years, season = month)
    x.september <- x.september0
    x.september$Dates$start <- gsub('-09-', '-01-', x.september$Dates$start)
    x.september$Dates$end <- gsub('-09-', '-01-', x.september$Dates$end)

    y.september0 <- subsetGrid(y, years = training.years, season = month)
    y.september <- y.september0
    y.september$Dates$start <- gsub('-09-', '-01-', y.september$Dates$start)
    y.september$Dates$end <- gsub('-09-', '-01-', y.september$Dates$end)

    y.newdata.september0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.september <- y.newdata.september0
    y.newdata.september$Dates$start <- gsub('-09-', '-01-', y.newdata.september$Dates$start)
    y.newdata.september$Dates$end <- gsub('-09-', '-01-', y.newdata.september$Dates$end)


    # October
    month <- 10

    x.october <- subsetGrid(x, years = training.years, season = 12)
    x.october0 <- subsetGrid(x, years = training.years, season = month)
    x.october$Data <- x.october0$Data

    y.october <- subsetGrid(y, years = training.years, season = 12)
    y.october0 <- subsetGrid(y, years = training.years, season = month)
    y.october$Data <- y.october0$Data

    y.newdata.october <- subsetGrid(y.newdata, years = correction.years, season = 12)
    y.newdata.october0 <- subsetGrid(y.newdata, years = correction.years, season = month)
    y.newdata.october$Data <- y.newdata.october0$Data


    # November
    x.november <- subsetGrid(x, years = training.years, season = 11)
    y.november <- subsetGrid(y, years = training.years, season = 11)
    y.newdata.november <- subsetGrid(y.newdata, years = correction.years, season = 11)


    # Aggregate
    x.aggregate <- bindGrid(x.september, x.october, x.november, dimension = 'time')
    y.aggregate <- bindGrid(y.september, y.october, y.november, dimension = 'time')
    y.newdata.aggregate <- bindGrid(y.newdata.september, y.newdata.october, y.newdata.november, dimension = 'time')


    # Correction
    message('---------- Correcting SON season ----------')

    if(method == 'pqm'){
      SON.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       window = window,
                                       fitdistr.args = list(densfun = "gamma"),
                                       precipitation = precipitation,
                                       method = method,
                                       wet.threshold = threshold)

    } else if(method == 'isimip'){
      message(stringr::str_c('Correction type: ', isimip.type))
      SON.corrected0 <- isimip(x = x.aggregate,
                               y = y.aggregate,
                               newdata = y.newdata.aggregate,
                               type = isimip.type,
                               threshold = threshold)

    } else if(method == 'scaling' | method == 'dqm'){
      message(stringr::str_c('Correction type: ', scaling.type))
      SON.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       method = method,
                                       scaling.type = scaling.type)

    } else(
      SON.corrected0 <- biasCorrection(x = x.aggregate,
                                       y = y.aggregate,
                                       newdata = y.newdata.aggregate,
                                       window = window,
                                       precipitation = precipitation

    ,
                                       method = method,
                                       wet.threshold = threshold)
    )


    # Union grid
    september.corrected <- subsetGrid(SON.corrected0, season = 1)
    september.corrected$Dates$start <- as.Date(y.newdata.september0$Dates$start)
    september.corrected$Dates$end <- as.Date(y.newdata.september0$Dates$end)

    october.corrected <- subsetGrid(SON.corrected0, season = 12)
    october.corrected$Dates$start <- as.Date(y.newdata.october0$Dates$start)
    october.corrected$Dates$end <- as.Date(y.newdata.october0$Dates$end)

    november.corrected <- subsetGrid(SON.corrected0, season = 11)
    november.corrected$Dates$start <- as.Date(november.corrected$Dates$start)
    november.corrected$Dates$end <- as.Date(november.corrected$Dates$end)

    SON.corrected <- bindGrid(september.corrected, october.corrected, november.corrected, dimension = 'time')
    SON.corrected$Dates$start <- gsub('-12-', '-10-', SON.corrected$Dates$start)
    SON.corrected$Dates$end <- gsub('-12-', '-10-', SON.corrected$Dates$end)

    SON.corrected.no.NA <- suppressMessages(filterNA(SON.corrected))

    message('---------- Done ----------')




    # Union of corrected grids ----

    message('---------- Union grids ----------')

    corrected.grid <- bindGrid(DJF.corrected.no.NA, MAM.corrected.no.NA,
                               JJA.corrected.no.NA, SON.corrected.no.NA, dimension = 'time')
    corrected.grid$Dates$start <- as.Date(corrected.grid$Dates$start)
    corrected.grid$Dates$end <- as.Date(corrected.grid$Dates$end)

    message('---------- Done ----------')

    # Return corrected grid
    return(corrected.grid)

    # end ---
}
