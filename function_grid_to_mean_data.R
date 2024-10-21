grid_to_mean_data <- function(SpatialGridDataFrame){

  grid <- raster::stack(grid2sp(SpatialGridDataFrame))

  db1 <- raster::as.data.frame(grid, xy=FALSE)
  db2 <- colMeans(db1, na.rm = TRUE)
  db3 <- data.frame(date=names(db2), value=as.numeric(db2))

  db3$date <- gsub('X', '', db3$date)
  db3$date <- gsub('\\.', '-', db3$date)
  db3$date <- gsub('--03', '', db3$date)
  db3$date <- gsub('--04', '', db3$date)
  db3$date <- gsub('-12-00-00-GMT', '', db3$date)
  db3$date <- gsub('-00-00-00-GMT', '', db3$date)

  cleaned_db <- db3

  return(cleaned_db)

}
