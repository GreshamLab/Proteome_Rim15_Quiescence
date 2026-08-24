# This file contains the data analysis and visualization with Dan

#Tables
master_tidy_protein.csv - is a table that composted of 110,531 rows and 10 columns,
Columns 1-6 are informations for sample annotation, and column 7 - 10 are values from different levels of measurments.

Col 1- Protein IDs
Col 2- Gene names that encode the protein
Col 3- the genotype background of where the protein is measured: WT vs RIM15KO
Col 4- the nutritional conditions nutrient that the cells were cultures: C - carbon limitation, P - phosphorous limitation
Col 5- samplingTime: T00, T06, T16, T30
Col 6- replicate: Rep1, Rep2, Rep3
Col 7- raw ratio
Col 8- median normalized ratio
Col 9- raw intensity, added 1 as a sudocounts 
Col 10- normalized intensity based on Anders & Huber's method, based on raw intensity +1. 

#Plots