library(CluMSID)
library(xcms)
library(Spectra)
library(MetaboCoreUtils)
library(Rdisop)

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


c_mz <- 223.07243   
c_rt <- 0.96*60

ms2sub <- getSpectrum(ms2spectras, "precursor", c_mz, mz.tol = 0.01) #(5*mz)/1e6
ms2sub <- getSpectrum(ms2sub, "rt", c_rt, rt.tol = 10)
#ms2sub <- getSpectrum(ms2sub, "annotation", "compound_8_urine_pos_rep3.mzML")

if(length(ms2sub) > 1){
  intensitats <- c()
  for(i in seq(ms2sub)){
    idx <- ((substring(gsub(".*\\.","", accessSpectrum(ms2sub[[i]])[,1]), 1, 1)>1) & 
      (accessSpectrum(ms2sub[[i]])[,1] > 60))
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

c_mz <- 162.0561                 
#c_rt <- 9.39*60
sp_ms2list <- filterPrecursorMz(object = sp_xdata, mz = c_mz + 0.01 * c(-1, 1))
sp_ms2list <- filterRt(sp_ms2list, rt = c_rt + 30 * c(-1, 1))
length(sp_ms2list)
unique(basename(dataOrigin(sp_ms2list)))

options(scipen = -5)
plotSpectra(sp_ms2list, #main = sps$name,
            labels = function(z) format(mz(z)[[1L]], digits = 6),
            labelSrt = -30, labelPos = 2, labelOffset = 0.1)
options(scipen = 5)

c_frag <- c()
for(i in 328:nrow(db)){
  c_frag <- c(c_frag, unlist(strsplit(db$fragments[i], "; ")))
}
c_frag <- unique(c_frag)
if(db$polarity[nrow(db)] == "POS"){
  c_pol <- "[M+H]+"
}else{
  c_pol <- "[M-H]-"
}
c_frag_mz <- c()
for(i in 1:length(c_frag)){
  c_frag_mz <- c(c_frag_mz, unlist(mass2mz(getMolecule(c_frag[i])$exactmass, c_pol)))
}
for(i in 1:length(c_frag_mz)){
  sp_ms2list <- filterPrecursorMz(object = sp_xdata, mz = c_frag_mz[i] + 0.01 * c(-1, 1))
  sp_ms2list <- filterRt(sp_ms2list, rt = c_rt + 30 * c(-1, 1))
  if(length(sp_ms2list) >0){
    print(c_frag[i])
    print(c_frag_mz[i])
    print(length(sp_ms2list))
    print(unique(basename(dataOrigin(sp_ms2list))))
    print("####################################################")
  }
}

#############################################################################

db = read.csv("database.csv")
dbx = db[db$name == "3-Hydroxy-DL-kynurenine (fragment)", ]
dbx[,c(1:5, 8, 9)]
i = 1
ms2clu_i <- extractMS2spectra(
  paste0("mzML/", dbx$path[i], "/", dbx$file[i], ".mzML"),
  min_peaks = 2,
  recalibrate_precursor = FALSE)
ms2clu_i <- unlist(ms2clu_i)
ms2clu_i <- getSpectrum(ms2clu_i, "precursor", 206.04588, mz.tol = 0.01)
if(length(ms2clu_i) > 1){
  rts <- c()
  for(j in 1:length(ms2clu_i)){
    rts <- c(rts, ms2clu_i[[j]]@rt)
  }
  ms2clu_i <- ms2clu_i[[closest(72.1, rts)]]
} else if(length(ms2clu_i) == 1){
  rts <- ms2clu_i@rt
}
specplot(ms2clu_i)