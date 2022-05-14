library(CompoundDb)
(cdb <- CompDb("CompDb_inhouse_0.sqlite", flags = RSQLite::SQLITE_RW))
sps <- Spectra(cdb)
deleteSpectra(cdb, sps$spectrum_id[which(sps$compound_id == "M0028" & sps$adduct == "[M+CHO2-CO2]-")])
