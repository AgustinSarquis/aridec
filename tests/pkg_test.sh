#!/bin/bash

rm *tar.gz
rm -r aridec.Rcheck
R CMD build ../Rpkg/
R CMD check *tar.gz
