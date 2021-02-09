table(duplicated(gsub(" \\(fragment.*", "", db$name)))

library(xcms)
library(Rdisop)
library(MetaboCoreUtils)
library(MsCoreUtils)
library(Spectra)
library(CompoundDb)

db <- read.csv("database.csv")
db$adduct <- gsub("H-H2O-CO]", "H-HCOOH]", db$adduct)
db$adduct <- gsub("H-H2O-NH3]", "H-NH3-H2O]", db$adduct)
db$adduct <- gsub("H-hexose]", "H-Hexose-H2O]", db$adduct)
db$adduct <- gsub("\\(H2O)3-C2H2", "\\H2O-H2O-C2H4O (McLafferty)", db$adduct)
db$adduct <- gsub("-C2H4O*", "-C2H4O (McLafferty)", db$adduct)

i <- nrow(db)

c_fml <- db$formula[i]
c_add <- db$adduct[i]
c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
if(is.na(c_mz)){
  mi <- paste0(
    substr(db$adduct[i], 1, 4),
    substr(db$adduct[i], nchar(db$adduct[i])-1, nchar(db$adduct[i]))
  )
  c_mz <- unlist(mass2mz(getMolecule(db$formula[i])$exactmass, mi))
  losses <- unlist(strsplit(substr(db$adduct[i], 6, nchar(db$adduct[i])-2), "-"))
  losses <- gsub("hexose", "C6H10O5", losses)
  for(j in 1:length(losses)){
    c_mz <- c_mz - getMolecule(losses[j])$exactmass
  }
}

xdata <- readMSData(files = paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                    mode = "onDisk")
chr <- chromatogram(xdata, mz = c_mz + 0.01 * c(-1, 1))
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
c_rt <-   212.2464
dev.off()
plot(chr, xlim = c(c_rt - 50, c_rt + 50))
abline(v = c_rt)
sps <- xdata[[closest(c_rt, rtime(xdata)#, duplicates = "closest"
)]]
sps <- as.data.frame(sps)
sps <- sps[c(unlist(matchWithPpm(c_mz, sps$mz, ppm = 15)),
               unlist(matchWithPpm(c_mz + 1.003355, sps$mz, ppm = 15)),
               unlist(matchWithPpm((c_mz + 1.003355*2), sps$mz, ppm = 15))),]
write.table(sps, paste0("sirius/", db$abr[i], "_", db$polarity[i], "_FS.txt"), 
            row.names = FALSE, col.names = FALSE)


xdata <- Spectra(paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                 backend = MsBackendMzR())
ms2list <- filterPrecursorMz(object = xdata, mz = c_mz + 0.01 * c(-1, 1))
ms2list <- ms2list[closest(c_rt, rtime(ms2list))]
sps <- data.frame(
  mz = unlist(mz(ms2list)),
  int = unlist(intensity(ms2list))
)
write.table(sps, paste0("sirius/", db$abr[i], "_", db$polarity[i], "_MS2.txt"), 
            row.names = FALSE, col.names = FALSE)



plotSpectra(ms2list, #main = sps$name,
            labels = function(z) format(mz(z)[[1L]], digits = 6),
            labelSrt = -30, labelPos = 2, labelOffset = 0.1)

