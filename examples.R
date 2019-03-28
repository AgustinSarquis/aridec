library(aridec)

# Load all entries of the database as an R list. Make sure to change the path argument to the correct location in your computer
aridec=loadEntries(path='~/aridec/data/')

# How many entries are in the database?
length(aridec) 

# How can I see the first entry?
aridec[[1]]

# How can I see entries by name?
aridec[["Austin2006"]]

# List all entry names
names(aridec)

# How can I see the time-series data only?
aridec[[1]]$data

# How can I see the nested variables, e.g. latitude?
aridec[[1]]$SiteInfo$Coordinates$Latitude

# Plot data from one entry
# It runs out of colors and it shows only the first few references when there are too many variables.
plotEntry(entry=aridec[["Poca2014"]])
