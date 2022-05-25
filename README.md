# aridec: database of aridlands decomposition studies
This repository contains a database of litter decomposition studies for arid, semiarid and subhumid regions. 
All data with corresponding metadata are stored in the folder `data`.
The folder `Rpkg` contains source files of an R package that can be used to load database entries into an R session. 
The folder scripts includes useful code for different applications.
The folder test contains code to load new database entries and ensure files are suited to the database format.

To install, open R and run:
```R
install.packages("devtools")
devtools::install_github('AgustinSarquis/aridec/Rpkg')
```
