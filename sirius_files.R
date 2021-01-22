library(xcms)
library(Rdisop)
library(MetaboCoreUtils)
library(MsCoreUtils)
library(Spectra)
library(CompoundDb)

db <- read.csv("database.csv")
i <- nrow(db)

c_fml <- db$formula[i]
c_add <- db$adduct[i]

xdata <- readMSData(files = paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                    mode = "onDisk")
c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
chr <- chromatogram(xdata, mz = c_mz + 0.01 * c(-1, 1))
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
c_rt <- 54.3983
plot(chr)
abline(v = c_rt)
sps <- xdata[[closest(c_rt, rtime(xdata)#, duplicates = "closest"
)]]
sps <- as.data.frame(sps)
sps <- sps[c(unlist(matchWithPpm(c_mz, sps$mz, ppm = 15)),
               unlist(matchWithPpm(c_mz + 1.003355, sps$mz, ppm = 15)),
               unlist(matchWithPpm((c_mz + 1.003355*2), sps$mz, ppm = 15))),]
write.table(sps, paste0("sirius/", db$abr[i], "_FS.txt"), 
            row.names = FALSE, col.names = FALSE)


xdata <- Spectra(paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                 backend = MsBackendMzR())
ms2list <- filterPrecursorMz(object = xdata, mz = c_mz + 0.01 * c(-1, 1))
ms2list <- ms2list[closest(c_rt, rtime(ms2list))]
sps <- data.frame(
  mz = unlist(mz(ms2list)),
  int = unlist(intensity(ms2list))
)
write.table(sps, paste0("sirius/", db$abr[i], "_MS2.txt"), 
            row.names = FALSE, col.names = FALSE)
