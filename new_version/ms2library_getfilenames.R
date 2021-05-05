library(readxl)
library(Rdisop)
library(CompoundDb)
inj <- read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/standards_injections.xlsx")
cmps <- read_xlsx("C:/Users/garciaalom/Google Drive/laboratory/db_compounds.xlsx")


i.cmp <- "Neochlorogenic acid"
cmps$ID_cmp[cmps$compound == i.cmp]
i.inj <- inj[grep(cmps$ID_cmp[cmps$compound == i.cmp], inj$ID_cmp),]
i.inj <- i.inj[i.inj$chromatography == "plasma", ]
i.inj$path
i.inj$filename
unlist(mass2mz(getMolecule(cmps$formula[cmps$compound == i.cmp])$exactmass, c("[M-H]-", "[M+H]+")))

