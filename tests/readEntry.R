
install.packages("yaml") #Install this package only the first time you run this code

# Test first if you can read metadata
library(yaml)

metadataAustin2006=yaml.load_file(input = "~/Repos/aridec/data/Austin2006/metadata.yaml") #Change the directory path to your local path
metadataAustin2006 # If the previous line works without error you can now see the entry as an R list

# Test if you can read the data
dataAustin2006=read.csv("~/aridec/data/Austin2006/data.csv")
plot(dataAustin2006) #If you can see a plot, then R can read the data

