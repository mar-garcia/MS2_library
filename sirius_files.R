library(xcms)
library(Rdisop)
library(MetaboCoreUtils)
library(MsCoreUtils)
library(Spectra)
library(CompoundDb)

db <- read.csv("database.csv")
db$adduct <- gsub("H-H2O-CO]", "H-HCOOH]", db$adduct)
db$adduct <- gsub("H-hexose]", "H-Hexose-H2O]", db$adduct)
db$adduct <- gsub("\\(H2O)3-C2H2", "\\H2O-H2O-C2H4O (McLafferty)", db$adduct)
db$adduct <- gsub("-C2H4O*", "-C2H4O (McLafferty)", db$adduct)

i <- nrow(db)

c_fml <- db$formula[i]
c_add <- db$adduct[i]
if(c_add == "[M-H-CH3]-"){
  c_mz <- getMolecule(c_fml)$exactmass - 1.007276 - getMolecule("CH3")$exactmass
} else if(c_add == "[M+H-CH3]+"){
  c_mz <- getMolecule(c_fml)$exactmass + 1.007276 - getMolecule("CH3")$exactmass
} else if(c_add == "[M+H-CH3-CO]+"){
  c_mz <- getMolecule(c_fml)$exactmass + 1.007276 - getMolecule("CH3")$exactmass - 
    getMolecule("CO")$exactmass
} else if(c_add == "[M-H-CH4N]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H]-")) - 
    getMolecule("CH4N")$exactmass
} else if(c_add == "[M-H-H2O-CO2-C9H18]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-H2O-CO2]-")) - 
    getMolecule("C9H18")$exactmass
} else if(c_add == "[M-H-H2O-NH3-C4O2]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-H2O]-")) - 
    getMolecule("C4H3NO2")$exactmass
} else if(c_add == "[M-H-H2O-C2O2]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-H2O]-")) - 
    getMolecule("C2O2")$exactmass
} else if(c_add == "[M+H-H2O-C2H2O-C2HNO]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-H2O]+")) - 
    getMolecule("C2H2O")$exactmass - getMolecule("C2HNO")$exactmass
} else if(c_add == "[M+H-(H2O)2-CO]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-(H2O)2]+")) - 
    getMolecule("CO")$exactmass
} else if(c_add == "[M+H-(H2O)2-CO-C5H4N2O]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-(H2O)2]+")) - 
    getMolecule("CO")$exactmass - getMolecule("C5H4N2O")$exactmass
} else if(c_add == "[M-H-(H2O)2-CO2]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-(H2O)2]-")) - 
    getMolecule("CO2")$exactmass
} else if(c_add == "[M+H-(H2O)2-C2H2O]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-(H2O)2]+")) - 
    getMolecule("C2H2O")$exactmass
} else if(c_add == "[M-H-(H2O)2-C3H2O]-" | c_add == "[M-H-(H2O)2-CO-C2H2]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-(H2O)2]-")) - 
    getMolecule("C3H2O")$exactmass
} else if(c_add == "[M+H-(H2O)3-CO-C2H4]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-(H2O)3-CO]+")) - 
    getMolecule("C2H4")$exactmass
} else if(c_add == "[M+H-(H2O)4]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-(H2O)3]+")) - 
    getMolecule("H2O")$exactmass
} else if(c_add == "[M+H-(H2O)4-C3H6]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-(H2O)3]+")) - 
    getMolecule("H2O")$exactmass - getMolecule("C3H6")$exactmass
} else if(c_add == "[M-H-CO2-C2H3NO]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-CO2]-")) - 
    getMolecule("C2H3NO")$exactmass
} else if(c_add == "[M-H-CO2-C3H5NO]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-CO2]-")) - 
    getMolecule("C3H5NO")$exactmass
} else if(c_add == "[M-H-hexose-H2O]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-Hexose-H2O]-")) - 
    getMolecule("H2O")$exactmass
} else if(c_add == "[M+H-hexose-H2O]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-Hexose-H2O]+")) - 
    getMolecule("H2O")$exactmass
} else if(c_add == "[M+H-hexose-(H2O)2-CH2]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H-Hexose-H2O]+")) - 
    getMolecule("H2O")$exactmass*2 - getMolecule("CH2")$exactmass
} else if(c_add == "[M-H-(hexose)2]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H]-")) - 
    getMolecule("C6H10O5")$exactmass*2
} else if(c_add == "[M+H-C2H5NO3]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H]+")) - 
    getMolecule("C2H5NO3")$exactmass
} else if(c_add == "[M-H-C10H18O9]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H]-")) - 
    getMolecule("C10H18O9")$exactmass
} else if(c_add == "[M-H-C9H7N-CO]-"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M-H-CO]-")) - 
    getMolecule("C9H7N")$exactmass
} else if(c_add == "[M+H-H3NO-CH2N2]+"){
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, "[M+H]+")) - 
    getMolecule("H3NO")$exactmass - getMolecule("CH2N2")$exactmass
} else {
  c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
}

xdata <- readMSData(files = paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                    mode = "onDisk")
chr <- chromatogram(xdata, mz = c_mz + 0.01 * c(-1, 1))
chromPeaks(findChromPeaks(chr, param = CentWaveParam(peakwidth = c(2, 20))))
c_rt <-  79.7955   
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
            labels = function(z) format(mz(z)[[1L]], digits = 4),
            labelSrt = -30, labelPos = 2, labelOffset = 0.1)

