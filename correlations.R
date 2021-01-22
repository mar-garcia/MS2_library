library(Spectra)

db <- read.csv("database.csv")
spd <- DataFrame(
  msLevel = rep(2L, nrow(db)),
  polarity = rep(1L, nrow(db)),
  id = db$abr,
  name = db$name)

l_mz <- list()
l_intensity <- list()
## Assign m/z and intensity values.
for(i in 1:nrow(db)){
  c_fml <- db$formula[i]
  c_add <- db$adduct[i]
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
  c_rt <- db$RT[i]
  xdata <- Spectra(paste0("mzML/", db$file[i], ".mzML"), 
                   backend = MsBackendMzR())
  ms2list <- filterPrecursorMz(object = xdata, mz = c_mz + 0.01 * c(-1, 1))
  ms2list <- ms2list[closest(c_rt, rtime(ms2list))]
  sps <- data.frame(
    mz = unlist(mz(ms2list)),
    int = unlist(intensity(ms2list))
  )
  c_frag <- c(db$formula[i], unlist(strsplit(db$fragments[i], "; ")))
  c_frag_mz <- c()
  for(j in 1:length(c_frag)){
    c_frag_mz <- c(c_frag_mz, unlist(mass2mz(getMolecule(c_frag[j])$exactmass, "[M-H]-")))
  }
  idx <- unlist(matchWithPpm(c_frag_mz, sps$mz, ppm = 10))
  l_mz[[i]] <- sps$mz[idx][order(sps$mz[idx])]
  l_intensity[[i]] <- sps$int[idx][order(sps$mz[idx])]
}
spd$mz <- l_mz
spd$intensity <- l_intensity
spd <- Spectra(spd)
rm(l_mz, l_intensity, j, idx, sps, xdata, ms2list,
   c_add, c_fml, c_frag, c_frag_mz, c_mz, c_rt)

####################################

library(CluMSID)
i <- 1
c_fml <- db$formula[i]
c_add <- db$adduct[i]
c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
c_rt <- db$RT[i]
c_frag <- c(db$formula[i], unlist(strsplit(db$fragments[i], "; ")))
c_frag_mz <- c()
for(j in 1:length(c_frag)){
  c_frag_mz <- c(c_frag_mz, unlist(mass2mz(getMolecule(c_frag[j])$exactmass, "[M-H]-")))
}
ms2clu <- extractMS2spectra(paste0("mzML/", db$file[i], ".mzML"),
                  min_peaks = 2,
                  recalibrate_precursor = FALSE)
ms2clu <- unlist(ms2clu)
ms2clu <- getSpectrum(ms2clu, "precursor", c_mz, mz.tol = 0.01)
ms2clu <- getSpectrum(ms2clu, "rt", c_rt, rt.tol = 2)
ms2clu@id <- db$name[i]
ms2clu@annotation <- db$adduct[i]
idx <- unlist(matchWithPpm(c_frag_mz, ms2clu@spectrum[,1], ppm = 10))
ms2clu@spectrum <- matrix(ms2clu@spectrum[idx,], ncol = 2)

for(i in 2:nrow(db)){
  c_fml <- db$formula[i]
  c_add <- db$adduct[i]
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
  c_rt <- db$RT[i]
  c_frag <- c(db$formula[i], unlist(strsplit(db$fragments[i], "; ")))
  c_frag_mz <- c()
  for(j in 1:length(c_frag)){
    c_frag_mz <- c(c_frag_mz, unlist(mass2mz(getMolecule(c_frag[j])$exactmass, "[M-H]-")))
  }
  ms2clu_i <- extractMS2spectra(paste0("mzML/", db$file[i], ".mzML"),
                                min_peaks = 2,
                                recalibrate_precursor = FALSE)
  ms2clu_i <- unlist(ms2clu_i)
  ms2clu_i <- getSpectrum(ms2clu_i, "precursor", c_mz, mz.tol = 0.01)
  rts <- c()
  for(j in 1:length(ms2clu_i)){
    rts <- c(rts, ms2clu_i[[j]]@rt)
  }
  ms2clu_i <- ms2clu_i[[closest(c_rt, rts)]]
  ms2clu_i@id <- db$name[i]
  ms2clu_i@annotation <- db$adduct[i]
  idx <- unlist(matchWithPpm(c_frag_mz, ms2clu_i@spectrum[,1], ppm = 10))
  ms2clu_i@spectrum <- matrix(ms2clu_i@spectrum[idx,], ncol = 2)
  
  ms2clu <- c(ms2clu, ms2clu_i)
}
rm(i, j, ms2clu_i, idx, c_add, c_fml, c_frag, c_frag_mz, c_mz, c_rt, rts)

####################################

x_mz <- c(87.0089, 111.0090, 129.0200, 173.0097)
x_int <- c(38801.211, 7319.501, 385760.406, 10595.318)
write.table(cbind(x_mz, x_int), "x_MS2.txt", row.names = F, col.names = F)

x_spd <- DataFrame(
  msLevel = 2L,
  polarity = 1L,
  id = "x",
  name = "x")
x_spd$mz <- list(x_mz)
x_spd$intensity <- list(x_int)
x_spd <- Spectra(x_spd)
c_spd <- c(x_spd, spd)

ms2clu_i@id <- "x"
ms2clu_i@annotation <- "x"
ms2clu_i@precursor <- 1000
ms2clu_i@rt <- 500
ms2clu_i@polarity <- "x"
ms2clu_i@spectrum <- cbind(x_mz, x_int)

cbind(c_spd$name[order(compareSpectra(c_spd, ppm = 50)[,1], decreasing = T)], 
      compareSpectra(c_spd, ppm = 50)[,1][order(compareSpectra(c_spd, ppm = 50)[,1], 
                                              decreasing = T)]
)[-1,]
tmp <- matrix((getSimilarities(ms2clu_i, ms2clu)[order(
  getSimilarities(ms2clu_i, ms2clu), decreasing = T)]), ncol = 1)
rownames(tmp) <- names((getSimilarities(ms2clu_i, ms2clu)[order(
  getSimilarities(ms2clu_i, ms2clu), decreasing = T)]))
tmp
