## ---
## title: "Example to fit the function 'multiSidbFit'"
## author: " "
## date: "03/15/2021"
## ---

## imports
packs <- c('SoilR', 'sidb', 'yaml', 'parallel','R.utils', 'aridec')
sapply(packs, require, character.only = TRUE) ## be sure all the outputs are TRUE

## Function uploading 
userLocation <- '~/aridec/scripts/multiSidbFit/multiSidbFit' #Change location of the sidb folder
source.R <- file.path(userLocation)
sourceDirectory(source.R, modifiedOnly = TRUE, verbose = TRUE)

## Data-entry loading 
path <- "~/aridec/data/" ## set a correct file path
load_entries <- loadEntries(path)
db <- load_entries[["Arriaga2007"]]
str(db)



## Example 1
twopsFit_model <- multiSidbFit(db,
                       model = 'twopsFit',
                       initialCarbon = 100,
                       ts.col = 10:15,
                       inipars=c(1.11, 1.11, 0.23, 0.95))
str(twopsFit_model)
names(twopsFit_model)

## Example 2
threeppFit_model <- multiSidbFit(db,
                    model = 'threeppFit',
                    ts.col = 10:15,
                    inipars=c(0.05, 0.01, 0.001, 0.1, 0.1))
str(threeppFit_model)

## Example 3
twopfFit_model <- multiSidbFit(db,
                               model = 'twopfFit',
                               initialCarbon = 100,
                               ts.col = 10:15,
                               inipars=c(1, 0.5, 0.5, 0.5, 0.3))