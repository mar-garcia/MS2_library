library(Spectra)
library(readxl)
cmps <- read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/db_compounds.xlsx")
norm_int <- function(x, ...) {
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
}
label_fun <- function(x) {
  ints <- unlist(intensity(x))
  mzs <- format(unlist(mz(x)), digits = 6)
  mzs[ints < 10] <- ""
  mzs
}

load("~/MS2library_POS.RData")

spd <- DataFrame(
  msLevel = c(2L),
  polarity = c(1L),
  id = c("X"),
  name = c("X"))

## Assign m/z and intensity values.
spd$mz <- list(
  c(179.6892,179.7395,179.7601,179.8158,179.883,179.9508,215.0559,232.0847,
    233.0676,250.0944,251.0785,261.2006,271.5399,283.8734,341.1587,341.9513,
    354.6637,394.1632,403.174,415.1726,433.1838,464.3528,475.4868,482.7329,
    518.2347,530.2349
  ))
spd$intensity <- list(
  c(7850,8711,9441,8894,8440,7695,7357,9422,9744,268092,32053,6605,8805,21350,
    9358,20591,10755,12561,41131,48849,51010,20707,7644,7686,21265,14690
  ))
spd$precursorMz <- 550.2624
sps <- Spectra(spd, backend = MsBackendDataFrame())

# Normalize the intensities
sps <- addProcessing(sps, norm_int)
sps.ms2 <- addProcessing(sps.ms2, norm_int)

# Pairwise comparison of all spectra
#cormat <- compareSpectra(sps.ms2, ppm = 20)
#hm <- pheatmap::pheatmap(cormat)

# Compare all MS2 agains "i"
res <- compareSpectra(sps.ms2, sps, ppm = 20)
# Identify the best matching pair
idx <- order(-res)[1]

plotSpectraMirror(sps, sps.ms2[idx], tolerance = 0.001, 
                  labels = label_fun, labelPos = 2, labelOffset = 0.2,
                  labelSrt = -30, main = cmps$compound[cmps$ID_cmp == sps.ms2$name[idx]])
grid()




# Neutral losses ---------------------------------------------
sps.ms2_nl <- applyProcessing(sps.ms2)
mz(sps.ms2_nl@backend) <- mz(sps.ms2_nl) - precursorMz(sps.ms2_nl)

sps_nl <- applyProcessing(sps)
mz(sps_nl@backend) <- mz(sps_nl) - precursorMz(sps_nl)

# create similarity matrix
#similarity <- compareSpectra(sps.ms2_nl, tolerance = 0.005)
#colnames(similarity) <-rownames(similarity) <- sps.ms2_nl$name
#similarity_nl <- compareSpectra(sps.ms2_nl, tolerance = 0.005)
#colnames(similarity_nl) <- rownames(similarity_nl) <- sps.ms2_nl$name

# Plot the heatmap
#paletteLength <- 50
#myColor <- colorRampPalette(c("blue", "red"))(paletteLength)
#myBreaks <- c(seq(0,1, by = 1/50))
#pheatmap::pheatmap(similarity, color=myColor, breaks=myBreaks, main = "Dotproduct")
#pheatmap::pheatmap(similarity_nl, color=myColor, breaks=myBreaks, main = "NL-Dotproduct ver 2")

#plotSpectraMirror(sps.ms2_nl[1], sps.ms2_nl[2])


res <- compareSpectra(sps.ms2_nl, sps_nl, tolerance = 0.005)
idx <- order(-res)[1]

plotSpectraMirror(sps_nl, sps.ms2_nl[idx], tolerance = 0.01, 
                  labels = label_fun, labelPos = 2, labelOffset = 0.2,
                  labelSrt = -30, main = cmps$compound[cmps$ID_cmp == sps.ms2$name[idx]])
