pkgname <- "aridec"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('aridec')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MAP")
### * MAP

flush(stderr()); flush(stdout())

### Name: MAP
### Title: Creates a data frame with mean annual precipitation values of
###   the sites
### Aliases: MAP

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
MAP=MAP(database=aridec)



cleanEx()
nameEx("MAT")
### * MAT

flush(stderr()); flush(stdout())

### Name: MAT
### Title: Creates a data frame with mean annual temperature values of the
###   sites
### Aliases: MAT

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
MAT=MAT(database=aridec)



cleanEx()
nameEx("biome")
### * biome

flush(stderr()); flush(stdout())

### Name: biome
### Title: Creates a data frame with the ecosystem type of the sites
### Aliases: biome

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
biome=biome(database=aridec)



cleanEx()
nameEx("carbon")
### * carbon

flush(stderr()); flush(stdout())

### Name: carbon
### Title: Creates a data frame with the carbon content in litter samples
### Aliases: carbon

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
C=carbon(database=aridec)



cleanEx()
nameEx("coordinates")
### * coordinates

flush(stderr()); flush(stdout())

### Name: coordinates
### Title: Creates a data frame with the coordinates of the sites
### Aliases: coordinates

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
coor=coordinates(database=aridec)



cleanEx()
nameEx("countries")
### * countries

flush(stderr()); flush(stdout())

### Name: countries
### Title: Creates a data frame with the countries of the sites
### Aliases: countries

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
countries=countries(database=aridec)



cleanEx()
nameEx("elevation")
### * elevation

flush(stderr()); flush(stdout())

### Name: elevation
### Title: Creates a data frame with elevation values of the sites
### Aliases: elevation

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
elevation=elevation(database=aridec)



cleanEx()
nameEx("lignin")
### * lignin

flush(stderr()); flush(stdout())

### Name: lignin
### Title: Creates a data frame with the lignin content in litter samples
### Aliases: lignin

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
lignin=lignin(database=aridec)



cleanEx()
nameEx("loadEntries")
### * loadEntries

flush(stderr()); flush(stdout())

### Name: loadEntries
### Title: Load all entries of the aridec dataset
### Aliases: loadEntries

### ** Examples

aridec=loadEntries()



cleanEx()
nameEx("material")
### * material

flush(stderr()); flush(stdout())

### Name: material
### Title: Creates a data frame with the list of litter samples' plant
###   parts
### Aliases: material

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
material=material(database=aridec)



cleanEx()
nameEx("nitrogen")
### * nitrogen

flush(stderr()); flush(stdout())

### Name: nitrogen
### Title: Creates a data frame with the nitrogen content of litter samples
### Aliases: nitrogen

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
N=nitrogen(database=aridec)



cleanEx()
nameEx("onepFit")
### * onepFit

flush(stderr()); flush(stdout())

### Name: onepFit
### Title: Fits a one pool model to a time-series
### Aliases: onepFit

### ** Examples

aridec<-loadEntries(path='~/Documents/GitHub/aridec/data/')
entry=aridec[[20]]
a=onepFit(timeSeries = entry$timeSeries[,1:2],
initialCarbon=100)



cleanEx()
nameEx("plotEntry")
### * plotEntry

flush(stderr()); flush(stdout())

### Name: plotEntry
### Title: Plot individual entries of the aridec dataset
### Aliases: plotEntry

### ** Examples

aridec<-loadEntries(path='~/Documents/GitHub/aridec/data/')
plotEntry(entry=aridec[["Adair2017"]])



cleanEx()
nameEx("readEntry")
### * readEntry

flush(stderr()); flush(stdout())

### Name: readEntry
### Title: Read single entry of the aridec database
### Aliases: readEntry

### ** Examples

Adair2017=readEntry(path="~/Documents/GitHub/aridec/data/", entryName="Adair2017")



cleanEx()
nameEx("soilorder")
### * soilorder

flush(stderr()); flush(stdout())

### Name: soilorder
### Title: Creates a data frame with soil orders of the sites
### Aliases: soilorder

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
soilorder=soilorder(database=aridec)



cleanEx()
nameEx("species")
### * species

flush(stderr()); flush(stdout())

### Name: species
### Title: Creates a data frame with the species list of litter samples
### Aliases: species

### ** Examples

aridec<-loadEntries(path='~/Repos/aridec/data/')
species=species(database=aridec)



cleanEx()
nameEx("threeppFit")
### * threeppFit

flush(stderr()); flush(stdout())

### Name: threeppFit
### Title: Fits a three pool model with parallel structure to a time series
### Aliases: threeppFit

### ** Examples




cleanEx()
nameEx("threepsFit")
### * threepsFit

flush(stderr()); flush(stdout())

### Name: threepsFit
### Title: Fits a three pool model with series structure to a time series
### Aliases: threepsFit

### ** Examples




cleanEx()
nameEx("twopfFit")
### * twopfFit

flush(stderr()); flush(stdout())

### Name: twopfFit
### Title: Fits a two pool model with feedback structure to a time series
### Aliases: twopfFit

### ** Examples




cleanEx()
nameEx("twoppFit")
### * twoppFit

flush(stderr()); flush(stdout())

### Name: twoppFit
### Title: Fits a two pool model with parallel structure to a time series
### Aliases: twoppFit

### ** Examples




cleanEx()
nameEx("twopsFit")
### * twopsFit

flush(stderr()); flush(stdout())

### Name: twopsFit
### Title: Fits a two pool model with series structure to a time series
### Aliases: twopsFit

### ** Examples




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
