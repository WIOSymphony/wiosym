# nit is the number of iterations through which the raster will "expand", i.e.
# the moving window will move. You will need more iterations if you want to
# extrapolate to cells distant from the existing data.
# The win_size is the radius of the circle around a pixel from which
# the mean values are taken in the moving window. The best choice is the resolution
# in map units of your target raster. But a larger size will increase speed.


resint <- function(source, grid, nit = 100, win_size = 1000, method = "cubic"){

  source <- terra::project(source, grid, method = method)

  for(i in 1:nit){
    
    w <- terra::focalMat(source, win_size, "circle")
    w[w > 0] <- 1
    f <- terra::focal(source, w = w, fun = "mean", na.policy = "only", na.rm=TRUE)
    source[is.na(source)] <- f[is.na(source)]
    source[is.na(grid)] <- NA
    perc <- (i/nit)*100
    cat("\r", perc, "% complete")

  }

  return(source)
 
}