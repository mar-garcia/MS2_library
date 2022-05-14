library(CompoundDb)
(cdb <- CompDb("CompDb_inhouse_0.sqlite", flags = RSQLite::SQLITE_RW))
sps <- Spectra(cdb)
deleteSpectra(cdb, sps$spectrum_id[which(sps$compound_id == "M0004" & sps$adduct == "[M-H]-")])

