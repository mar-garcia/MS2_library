library(CluMSID)
library(xcms)

mzXMLfiles <- list.files("tmp/")
spectras <- lapply(paste0("tmp/", mzXMLfiles), 
                   function(x){
                     extractMS2spectra(x,
                                       min_peaks = 2,
                                       recalibrate_precursor = FALSE)
                   })

ms2list <- unlist(spectras)

muestra <- NA
for(i in 1:length(spectras)){
  muestra <- c(muestra,
               rep(mzXMLfiles[[i]],
                   length(spectras[[i]])))
}
muestra <- muestra[!is.na(muestra)]
for(i in 1:length(muestra)){
  slot(ms2list[[i]], "annotation") <- gsub(".*\\/", "", muestra[i])
}


c_mz <- 191.0197
c_rt <- 1.2*60

ms2sub <- getSpectrum(ms2list, "precursor", c_mz, mz.tol = 0.001) #(5*mz)/1e6
ms2sub <- getSpectrum(ms2sub, "rt", c_rt, rt.tol = 10)

if(length(ms2sub) > 1){
  intensitats <- c()
  for(i in seq(ms2sub)){
    idx <- substring(gsub(".*\\.","", accessSpectrum(ms2sub[[i]])[,1]), 1, 1)>1
    int.noise <- accessSpectrum(ms2sub[[i]])[idx,2][which.max(accessSpectrum(ms2sub[[i]])[idx,2])]
    int.good <- accessSpectrum(ms2sub[[i]])[-idx,2][which.max(accessSpectrum(ms2sub[[i]])[-idx,2])]
    intensitats <- c(intensitats, int.good / int.noise)
  }
}

dev.off()
par(mfrow=c(1,2))
if(length(ms2sub) > 30){
  for(i in (length(ms2sub)-30):length(ms2sub)){
    j <- order(intensitats)[i]
    raw_data <- readMSData(files = paste0("tmp/", ms2sub[[j]]@annotation), 
                           mode = "onDisk")
    chr <- chromatogram(raw_data, 
                        mz = c_mz + 0.01 * c(-1, 1), 
                        rt = c_rt + 20 * c(-1, 1))
    plot(chr, xlim = c_rt + 20 * c(-1, 1))
    abline(v=ms2sub[[j]]@rt)
    specplot(ms2sub[[j]])
  }
} else if(length(ms2sub) > 1 & length(ms2sub) <= 30){
  for(i in 1:length(ms2sub)){
    j <- order(intensitats)[i]
    raw_data <- readMSData(files = paste0("tmp/", ms2sub[[j]]@annotation), 
                           mode = "onDisk")
    chr <- chromatogram(raw_data, 
                        mz = c_mz + 0.01 * c(-1, 1), 
                        rt = c_rt + 50 * c(-1, 1))
    plot(chr, xlim = c_rt + 50 * c(-1, 1))
    abline(v=ms2sub[[j]]@rt)
    specplot(ms2sub[[j]])
  }
} else if(length(ms2sub) == 1){
  raw_data <- readMSData(files = paste0("tmp/", ms2sub@annotation), 
                         mode = "onDisk")
  chr <- chromatogram(raw_data, 
                      mz = c_mz + 0.01 * c(-1, 1), 
                      rt = c_rt + 20 * c(-1, 1))
  plot(chr, xlim = c_rt + 20 * c(-1, 1))
  abline(v=ms2sub@rt)
  specplot(ms2sub)
}