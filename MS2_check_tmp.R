library(CluMSID)
library(xcms)
library(Spectra)

mzXMLfiles <- list.files("tmp/")
spectras <- lapply(paste0("tmp/", mzXMLfiles), 
                   function(x){
                     extractMS2spectra(x,
                                       min_peaks = 2,
                                       recalibrate_precursor = FALSE)
                   })

ms2spectras <- unlist(spectras)

muestra <- NA
for(i in 1:length(spectras)){
  muestra <- c(muestra,
               rep(mzXMLfiles[[i]],
                   length(spectras[[i]])))
}
muestra <- muestra[!is.na(muestra)]
for(i in 1:length(muestra)){
  slot(ms2spectras[[i]], "annotation") <- gsub(".*\\/", "", muestra[i])
}


c_mz <- 181.0718
c_rt <- 0.82*60

ms2sub <- getSpectrum(ms2spectras, "precursor", c_mz, mz.tol = 0.1) #(5*mz)/1e6
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


#############################################

i <- 1
sp_xdata <- Spectra(paste0("tmp/", mzXMLfiles[i]), 
                 backend = MsBackendMzR())
for(i in 2:length(mzXMLfiles)){
  sp_xdata <- c(sp_xdata, Spectra(paste0("tmp/", mzXMLfiles[i]), 
                            backend = MsBackendMzR()))
}

c_mz <- 101.0245
c_rt <- 0.82*60
sp_ms2list <- filterPrecursorMz(object = sp_xdata, mz = c_mz + 0.01 * c(-1, 1))
sp_ms2list <- filterRt(sp_ms2list, rt = c_rt + 10 * c(-1, 1))
length(sp_ms2list)