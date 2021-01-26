library(xcms)
library(Rdisop)
library(MetaboCoreUtils)
library(MsCoreUtils)
library(Spectra)
library(CompoundDb)

db <- read.csv("database.csv")
db$adduct <- gsub("M-H-hexose", "M-H-Hexose-H2O", db$adduct)
i <- nrow(db)

c_fml <- db$formula[i]
c_add <- db$adduct[i]
if(c_add == "[M-H-CH3]-"){
  c_mz <- getMolecule(c_fml)$exactmass - 1.007276 - getMolecule("CH3")$exactmass
} else if(c_add == "[M-H-(H2O)2-CO2]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-(H2O)2]-")) - 43.98982926
} else {
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
}

xdata <- readMSData(files = paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                    mode = "onDisk")
chr <- chromatogram(xdata, mz = c_mz + 0.01 * c(-1, 1))
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
c_rt <- 483.2010
dev.off()
plot(chr, xlim = c(c_rt - 30, c_rt + 30))
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
            labels = function(z) format(mz(z)[[1L]], digits = 4),
            labelSrt = -30, labelPos = 2, labelOffset = 0.1)
