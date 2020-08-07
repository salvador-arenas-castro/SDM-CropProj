
convertToGeoTIFF <- function(inputFolder,outputFolder){
  
  fl <- list.files(inputFolder, pattern= ".grd$", full.names = TRUE)
  
  for(f in fl){
    
    cat("Converting file:\n",f," ...\n",sep="")
    fn <- tools::file_path_sans_ext(basename(f))
    
    raster::stack(f) %>%
      raster::writeRaster(filename = paste(outputFolder,"/",fn,".tif",sep=""))
    
    cat("done!\n\n")
  }
}


getPercTrain <- function(Np, maxTrainPoints = 3000, defaultPerc = 80){
  
  Nall <- 2*Np # Assumes that sample is double the nr of presences
  
  if(Nall > maxTrainPoints){
    perc <- ceiling((maxTrainPoints / Nall)*100)
    return(perc)
    
  }else{
    return(defaultPerc)
  }
}
