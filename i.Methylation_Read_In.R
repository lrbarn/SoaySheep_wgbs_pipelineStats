#### i.Methylation_Read_In.R####
## This gets the samples list from the bismark files directory and then sets up the methylation database before filtering (ii) and then loocv (iii)

### NOTE
# make sure you are in the right working directory before starting!
setwd("/Users/bip23lrb/Documents/Sheffield/Projects/PracticeThings/MethylKit/")

#### Packages ####
library(tidyverse)
library(methylKit)

print("starting methylation import")

#### Data ####
sample.list <- c(list.files("Data/"))
# lists all the sample file names
# need the full path
full_path <- paste0("Data/", sample.list[1:length(sample.list)])

#convert the resulting character vector into a list for methRead
file.list <- as.list(full_path)

## making a list of ids
sampleIDs <- c(str_remove(sample.list, "\\..*$"))
id.list <- as.vector(sampleIDs)
id.list <- as.list(sampleIDs)

## treatment vector
treatments <- c(rep(0,length(sample.list)))
treatment.list <- as.vector(treatments)
treat.vec <- c(0,0,0)

### Importing the data ###
methylationDB <- methRead(file.list,
                          sample.id = id.list,
                          treatment = treat.vec,
                          context = "CpG",
                          dbtype = "tabix",
                          dbdir = "sampleDB",
                          pipeline = "bismarkCoverage",
                          #mincov = 10,
                          assembly = "sheep"
                          )

print("Methylation import completed")

#### Initial Filtering ####
## this is being kept to a minimum to preserve as many CpGs as possible
filtered.methylationDB <- filterByCoverage(methylationDB,
                                           #lo.count = 10,
                                           hi.perc = 99)

print("methylation filtering completed")

methylationCpGs <- methylKit::unite(filtered.methylationDB, destrand = TRUE)

print("methylation data united")

# goal is to get the methylationCpGs output into the right format to put into the elastic net

print("Reformatting methylation data")

# extracing the simple data out of the methylKit output
simple <- getData(methylationCpGs)
# extracting the sampleIDs that matches up in order to the simple df
Ids <- methylationCpGs@sample.ids

#### Creating a coordinates for each CpG
site <- paste(simple$chr, simple$start, sep = "_")

#### Creating a for loop to group together the output into useable data ####
## this loop is for proportion methylated per sample per site

# empty
columnLocation <- NA
IdLocation <- NA

# dataframe to put the output into
Prop_methylated <- data.frame(site = site)

#for loop to combine
for(i in seq(5, ncol(simple), 3)) { # bounces through each coverage column
  cov <- simple[,i] #list of coverage
  methCount <- simple[,i+1] # list of methylated reads
  methProp <- methCount/cov
  
  columnLocation <- ((i-2)/3) + 1
  IdLocation <- ((i-2)/3)
  
  Prop_methylated[columnLocation] <- methProp
  names(Prop_methylated)[columnLocation] <- Ids[IdLocation]
}

write_csv(Prop_methylated, "OutputData/Prop_methylated.csv")

print("Prop_methylated.csv written")
