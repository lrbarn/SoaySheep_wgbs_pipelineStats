#### Basis Script for Methylation Data Extraction ####

#===========NOTE==============
# there are places where the script will need to be adjusted depending on where the final data is stored
# this points of adjustment are demarked with        
#=============================

#### Packages ####
library(tidyverse)
library(methylKit)

#### Data and Directories ####
setwd() # this will depend on where the analysis is carried out

## for example doing it in my shared area of Stanage
# /users/bip23lrb/methylation_downstream
# can create a sensible file structure within this for where I want the data to be

sample.list <- c(list.files("Data/"))
#listing all the sample file names in the Data/ directory
#this is not the full file path, this will be adjusted below

# this will need to be adjusted depending on where the data's final location is             

file.list <- list(paste0("Data/", sample.list[1]),
                  paste0("Data/", sample.list[2]))
# this is clunky and will need to be adjusted again
# the goal here is to make a LIST of the file paths for every ".merged.bismark.cov.gz" file that will need to be read in
#       

id.list <- c(str_remove(sample.list, "\\..*$"))
# this removes the extension from the sample file names leaving just the sample id
# in this case its 000-000 which is our id and liverpool id
# could split again to just have one or the other if necessary
# this needs to be converted into a LIST FORMAT
#       
# an alternative is to just use list("id", "id", etc) but this is not going to work in the case of 800+ samples

treatment.vec <- c(rep(0, length(id.list)))

# Sample Annotation
# this is a dataframe of all the processing metadata about the samples eg. extraction group/ date, sequencing batch etc.
sampleAnnotation <- data.frame(batch_id =    ,
                               extraction_date =    )
## likely this will be in a dataframe that can just be imported but need to make sure the order matches up with the file list

#===========================================
#### METHYL KIT ####
#===========================================
####Importing the data #####
methylationDB <- methRead(file.list, # list of raw files
                          sample.id = id.list, # list of IDs
                          assembly = "Sheep", # label only
                          treatment = treatment.vec, # treatments, all 0
                          context = "CpG", # only looking at CpG sites
                          dbtype = "tabix", # save the output into a flat file database since it will be huge
                          dbdir = "methylationDB", # name of a new directory where the database will be stored
                          pipeline = "bismarkCoverage", # this specifies that the data is from a bismark format
                          mincov = 10 # not necessary at this point but can set to a minimun        
                          )

#### Descriptive stats #####
## methylation 
getMethylationStats(methylationDB[[1]], plot = TRUE, both.strands = FALSE)
# could change both.strands = TRUE since working on CpG data
# this is only done on a per sample basis when plot = TRUE
# not sure what happens when trying to call hundreds...

## Coverage
getCoverageStats(methylationDB[[2]], plot = TRUE, both.strands = FALSE)
# same notes as for methylation

#### Coverage Filtering #####
filtered.methylationDB <- filterByCoverage(methylationDB,
                                           lo.count = 10, # removes anything with a coverage less than 10x
                                           hi.perc = 99.9 # removes anything above the 99.9th percentile to account for PCR bias
                                           )
# can then do the sats again on the filtered data if you want
#     

#### Merging into the useable dataframe #####
methylationCpGs <- unite(filtered.methylationDB, destrand = TRUE)
# check what the destrand means for methylation data     
# this puts the data into the final format for some analysis
# it will also drop any CpG that is not present in all of the samples at this point, meaning there should be no gaps and each cpg should meet a coverage criteria

### More descriptive stats #####
getCorrelation(methylationCpGs, plot = FALSE)
# this will give the correlation between all samples
# the full data may be a bit large to do this for 
# plot = TRUE will make a nice graphic
# could do for just some samples      

clusterSamples(methylationCpGs, plot = TRUE)
# this creates a dendrogram of the sample clustering
# not sure how to manipulate this but could be a fun graphic?

### Batch Effects #####
# this is where the sample annotation dataframe is needed
assocComp(mBase = methylationCpGs,
          sampleAnnotation = sampleAnnotation)

#### Finding DMRs #####
# will find DMRs between the test and control groups
# not sure if this will work for individuals but can be a good place to start between sexes etc
DMRs <- calculateDiffMeth(methylationCpGs, mc.cores = 2)
# this takes a while to run....
# but can be split over multiple cores using mc.cores = 2 argument

## can then see where they are across the chromosomes or whether they are hypo/hyper methylated 
diff25p.hyper <- getMethylDiff(DMRs, difference =25, qvalue = 0.01, type = "hyper")
diff25p.hypo <- getMethylDiff(DMRs, difference = 25, qvalue = 0.01, type = "hypo")

#### Correcting Overdispersion ####
disperseDMR <- calculateDiffMeth(DMRs, overdispersion = "MN",
                                 test = "Chisq",
                                 mc.cores = 2)

#========================================
#========Prepping for Elastic Net =======
#========================================

#### Formatting the Data ####
#### Extracting the Methylkit data #####
methData <- getData(fullCpGs) 
# getting the table out of the united function from the methylkit part of the script

Ids <- fullCpGs@sample.ids
# get the sample IDs from the methylkit object
# do this here to make sure that the order of output and IDs matches up correctly

## creating coordinates for each CpG
# pastes the chromosome and the start postition of the base
# separated by "_" to make it possible to break down again if needed
site <- paste0(methData$chr, methData$start, sep = "_")

#### For loop to format into useful shape #####
## creates a dataframe of the proportion methylated per sample per site

#empty variables
columnLocation <- NA
IdLocation <- NA

#dataframe to put the output into
Prop_methylated <- data.frame(site = site)

# the loop
for (i in seq(5, ncol(simple), 3)) { #bounces through each coverage column
  cov <- methData[,i] #coverage for ith sample
  methCount <- methData[,i+i] # methylated reads for ith sample
  methProp <- methCount/cov
  
  columnLocation <- ((i-2)/3) + 1 # coordinates so the proportion goes into the correct column in the output df
  IdLocation <- ((i-2)/3) # finding the matching ID from the ID list to name the column the right thing
  
  #putting the data into the output dataframe
  Prop_methylated[columnLocation] <- methProp
  names(Prop_methylated)[columnLocation] <- Ids[IdLocation]
}

#### Writing into a CSV #####
write_csv(Prop_methylated, "OutputData/Prop_methylated.csv")