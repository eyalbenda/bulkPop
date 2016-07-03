#' @import parallel
initiateCluster<-function(ncores=7)
{
  tryCatch(stopCluster(cl),error=function(e)print("starting new cluster"))
  assign(x = "cl",makeCluster(getOption("cl.cores", ncores)),envir = .GlobalEnv)
}
