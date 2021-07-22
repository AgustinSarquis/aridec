# This script automatically generates the documentation of the aridec package

if(!require("devtools")) install.packages("devtools")

devtools::document("../Rpkg/")

