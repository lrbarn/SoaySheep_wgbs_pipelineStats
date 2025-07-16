#### Libraries ####
library(tidyverse)
library(methylKit)
library(glmnet)
# remember there will be some overlap in funciton names

#### Data ####
sample.list <- c(list.files("bismarkFiles/"))
# lists all the sample file names
# Generate the full paths for the first 7 files
full_paths <- paste0("bismarkFiles/", sample.list[1:7])

# Convert the resulting character vector into a list
file.list.alt <- as.list(full_paths)

# lists all the file paths for the samples
# needed for creating the methylation database

# taking the beginning of each file name as the sampleID
id.list <- c(str_remove(sample.list, "\\..*$"))
#making a list of sample IDs in the same order as the files
id.list.alt <- list("001-1", "104-108", "386-403", "483-502", "803-837", "821-855", "828-862")


#### Importing the data #####
methylationDB <- methRead(file.list.alt,
                          sample.id = id.list.alt,
                          assembly = "Sheep", # THIS NEEDS CHANGING TO CORRECT ASSEMBLY but just acts as a label
                          treatment = c(0,0,0,0,0,0,0),
                          context = "CpG",
                          dbtype = "tabix",
                          dbdir = "sampleDB",
                          pipeline = "bismarkCoverage",
                          mincov = 10
)


getMethylationStats(methylationDB[[2]], plot = FALSE, both.strands = FALSE)

filtered.methylationDB <- filterByCoverage(methylationDB,
                                           lo.count = 10,
                                           hi.perc = 99.9)

methylationCpGs <- unite(filtered.methylationDB, destrand = TRUE)

# goal is to get the methylationCpGs output into the right format to put into the elastic net


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

### now onto the elastic net
### BIG GAP TO FILL IN HERE ###

#### LOOCV CLOCK ####
## Data Prep
ID_list <- as.vector(unique(CpGmethylation$ID))
n <- length(ID_list)

nonZeroCpGs_names <- NA
variable_names <- c(Prop_methylated$site)

## admin
set.seed(123)

i <- NA
pb <- progress_bar$new(total = n)

training <- NA
test <- NA
predictors <- NA
response <- NA

## model set up
TRAIN_glmnet_CV <- NA
coef.lambda.min <- NA
run <- NA
testAge <- NA
testAgePRED <- NA
lambda.min <- NA
coef_min <- NA
nonZero <- NA
testID <- NA
nonZero_names_list <- list(NA)

#### The LOAOCV Loop ####
##     IN THEORY THIS WORKS BUT NEEDS A FOURTH SAMPLE SINCE THE MIN IS 3
for (i in 1:n) {
  #1 Setting up the training and test data (one id in test)
  training <- CpGmethylation %>% filter(ID != ID_list[i]) %>% dplyr::select(-ID)
  test <- CpGmethylation %>% filter(ID == ID_list[i]) %>% dplyr::select(-ID)
  
  #2 Setting up response and predictors
  predictors <- as.matrix(training %>% dplyr::select(-AgeY))
  response <- as.matrix(training %>% dplyr::select(AgeY))
  
  #3 Fitting the CV model
  TRAIN_glmnet_cv <- cv.glmnet(predictors, response, alpha = 0.5, type.measure = "mse") 
    ### WARNING MESSAGE HERE IS OK (currently only using 7 samples)
  #4 storing the lambda min for the run
  lambda.min[i] <- TRAIN_glmnet_cv$lambda.min
  
  #5 running again for the lambda min
  TRAIN_glmnet <- glmnet(predictors, response, alpha = 0.5, nlambda = 100)
  coef.lambda.min <- coef(TRAIN_glmnet)[, TRAIN_glmnet$lambda == lambda.min[i]]
  
  #6 Predicting the age(s) for the left out individual
  AgePrediction <- predict(TRAIN_glmnet, newx = (as.matrix(test %>% dplyr::select(-AgeY))), type = "response", s = lambda.min[i])
  
  #7 Collating the non zero coefficients
  coef_min <- coef(TRAIN_glmnet, s = lambda.min[i])
  nonZero[i] <- sum(coef_min[-1,] != 0)
  
  nonZero_names <- variable_names[coef_min[-1,] != 0]
  nonZero_names_list[[i]] <- nonZero_names
  
  #8 saving the outputs
  run[i] <- i
  testAge <- append(testAge, test$AgeY)
  testAgePRED <- append(testAgePRED, as.numeric(AgePrediction))
  testID <- c(testID, rep(ID_list[i], times = (length(AgePrediction))))
  
  # update the progress bar
  pb$tick()}

# collating output
LOAOCV_output <- data.frame(runNumber = run,
                            lambda.min = lambda.min,
                            nonZero = nonZero)
LOAOCV_predictions <- data.frame(testAge = testAge,
                                 testAgePRED = testAgePRED,
                                 #testID = testID,
                                 error = testAge - testAgePRED,
                                 sq_error = (testAge - testAgePRED)^2,
                                 ab_error = abs(testAge - testAgePRED),
                                 testID = testID) %>% 
  drop_na()

# pulling together the frequency of each coefficient
variable_freq <- as.data.frame(sort(table(unlist(nonZero_names_list)), decreasing = TRUE))
colnames(variable_freq) <- c("VarName", "Freq")


#### Plotting the results of LOAOCV####
# The individual age against predicted age based on the left out model
clockEstimates <- 
ggplot(LOAOCV_predictions, aes(x = testAge, y = testAgePRED)) +
  geom_point() + 
  labs(x = "Age (years)",
       y = "Predicted Age (Phenotype)",
       title = "Epigenetic Clock 1.0") +
  xlim(-0.5, 11) +
  ylim(-0.5, 11) +
  geom_abline(intercept = 0, slope = 1, colour = "pink", linetype = 2) +
  theme_bw()