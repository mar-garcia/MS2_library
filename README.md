# MS2 library
This repo contains a shiny app and additional codes for creating and exploring an in-house MS2 library.


## Main Steps
  
1. Transform the `.raw` DDA files to `.mzML` format with `Proteowizard`.   
2. Save the files in the folder `mzML` using the same path in which the original files are stored.      
3. Add in the file `database.csv` the info corresponding to the columnes: "name", "abr", "formula", "path", "file", "adduct".      
4. Apply the code `sirius_files.R` to check the RT (and insert this info in the file "database.csv") and to create the txt files of isotopic pattern and MS2 spectra.  
5. Using the software `Sirius` check which ions in the MS2 spectra match with the compound and add their formulas in the column "fragments" of the files "database.csv".  
6. Apply the code `MS2_library.Rmd` to update the inhouse database.  
7. Run the shiny app (`app.R`) to interactively use the inhouse database.  
  