library(plotly)

db <- read.csv("database.csv")
i <- 2

c_fml <- db$formula[i]
c_add <- db$adduct[i]

xdata <- readMSData(files = paste0("mzML/", db$file[i], ".mzML"), 
                    mode = "onDisk")
c_mz <- unlist(mass2mz(getMolecule(c_fml)$exactmass, c_add))
chr <- chromatogram(xdata, mz = c_mz + 0.01 * c(-1, 1))
c_rt <- db$RT[i]
c_chr <- data.frame(
  rt = rtime(chr[[1]])/60,
  int = intensity(chr[[1]])
)
c_chr[is.na(c_chr)] <- 0


xdata <- Spectra(paste0("mzML/", db$file[i], ".mzML"), 
                 backend = MsBackendMzR())
ms2list <- filterPrecursorMz(object = xdata, mz = c_mz + 0.01 * c(-1, 1))
rts <- rtime(ms2list)
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
sps$color <- "x"
sps$color[idx] <- "frag"
sps$color[which(sps$color == "x")[1]] <- "y"


plot_ly(data = c_chr, x= ~rt, y = ~int, mode = "lines") %>%
  add_annotations(x = rts/60, 
                  y = c_chr$int[closest(rts/60, c_chr$rt)], "x", showarrow = F) %>%
  add_annotations(x = rtime(ms2list)/60, 
                  y = c_chr$int[closest(rtime(ms2list)/60, c_chr$rt)], "MS2")
plot_ly() %>%
  add_segments(x = sps$mz, xend = sps$mz, 
               y = rep(0,nrow(sps)), yend = sps$int,
               color = sps$color) %>% 
  layout(showlegend = FALSE) %>%
  add_annotations(x = sps$mz[idx],
                  y = sps$int[idx],
                  round(sps$mz[idx], 4))
