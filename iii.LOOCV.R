#### iii.LOOCV.R ####
## for carrying out LOOCV with data read in by i.Methylation_Read_In and filtered by ii.Variance_Filtering
## follows on from ii.Variance_Filtering
print("=========================================")
print("Starting iii.LOOCV.R")
print("This version is LOanimalOCV with AgeY")
print("=========================================")

runName <- "pilot7"

print(paste0("Current run name: ", runName))

#### Packages ####
library(tidyverse)
library(methylKit)
library(glmnet)
library(progress)

#### Data ####
prop_methylated <- read_csv("data/Prop_methylated.csv")

#### Reformatting #####
## need columns as sites and rows as ids
#pivot long
prop_methylated_long <- prop_methylated %>%
  pivot_longer(
    cols = -site,            # Select all columns EXCEPT 'site'
    names_to = "sampleID",     # New column 'Sample' will store the original sample names
    values_to = "Methylation_Value" # New column 'Methylation_Value' will store the numeric values
  )
#pivot wide
CpGmethylation <- prop_methylated_long %>%
  pivot_wider(
    names_from = site,           # Values from the 'site' column become new column names
    values_from = Methylation_Value # Values from 'Methylation_Value' fill the new columns
  )



### Adding in Age 
## importing the metadata file
# metaData <- read_csv("~/Documents/Sheffield/Projects/DataWrangling/Reformatted_Data/metadata_all_REFORMAT.csv")
# metaData <- metaData %>% dplyr::select(sampleID, `CapYear-BirthYear`) %>% 
#   mutate(AgeY = `CapYear-BirthYear`) %>% 
#   dplyr::select(-`CapYear-BirthYear`)
# 
# metaData <- metaData %>% 
#   filter(sampleID %in% names(prop_methylated))
# 
# CpGmethylation <- left_join(CpGmethylation, metaData, by = join_by(sampleID))

##         THIS NEEDS A MAJOR FIX AND THEN TO MATCH UP WITH THE HPC VERSION        

print("Data formatting completed")

### now onto the elastic net

print("setting up the LOOCV and elastic net")

#### LOOCV CLOCK ####
## Data Prep
ID_list <- as.vector(unique(CpGmethylation$ID))
n <- length(ID_list)

nonZeroCpGs_names <- NA
variable_names <- c(prop_methylated$site)

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

print("LOOCV loop started")

#### The LOOCV Loop ####
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
  print(paste0(i, " of ", n))
  pb$tick()}

print("LOOCV loop completed")
print("output collation started")

# collating output
LOAOCV_output <- data.frame(runNumber = run,
                            lambda.min = lambda.min,
                            nonZero = nonZero,
                            testID = ID_list)
LOAOCV_predictions <- data.frame(testID = testID,
                                 testAgeY = testAge,
                                 testAgeYPRED = testAgePRED,
                                 error = testAge - testAgePRED,
                                 sq_error = (testAge - testAgePRED)^2,
                                 ab_error = abs(testAge - testAgePRED)) %>% 
  drop_na()

# pulling together the frequency of each coefficient
variable_freq <- as.data.frame(sort(table(unlist(nonZero_names_list)), decreasing = TRUE))
colnames(variable_freq) <- c("CpG", "Freq")

MeanSE <- mean(LOAOCV_predictions$sq_error)
MedianAbsoluteError <- median(LOAOCV_predictions$ab_error)
LOAOCV_correlation <- cor.test(LOAOCV_predictions$testAgeY, LOAOCV_predictions$testAgeYPRED)

print("Output collation done")
print("Creating Plots")

#### Plotting the results of LOAOCV####
# The individual age against predicted age based on the left out model
clockEstimates <- 
  ggplot(LOAOCV_predictions, aes(x = testAgeY, y = testAgeYPRED)) +
  geom_point() + 
  labs(x = "Age (years)",
       y = "Predicted Age (Phenotype)",
       title = "Epigenetic Clock 1.0") +
  xlim(-0.5, 11) +
  ylim(-0.5, 11) +
  geom_abline(intercept = 0, slope = 1, colour = "pink", linetype = 2) +
  theme_bw()

# number of non-zero CpG coefficients selected by models
numberSelectedCpGs <- 
  ggplot(LOAOCV_output, aes(x = nonZero)) +
  geom_bar() +
  labs(x = "Number of Selected CpGs") +
  theme_bw()

# Frequency each CpG is selected
CpGFrequency <-
  ggplot(variable_freq, aes(x = CpG, y = Freq)) +
  geom_col() +
  labs(x = "CpG",
       y = "Selection Frequency") +
  geom_text(data = variable_freq, aes(x = CpG, y = (Freq + 2), angle = 45),
            label = paste("n = ", variable_freq$Freq)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print("Plots Created")

#### pdf creation of plots ####
ggsave("clockEstimates.pdf", plot = clockEstimates, dpi = 300)
ggsave("numberSelectedCpGs.pdf", plot = numberSelectedCpGs, dpi = 300)
ggsave("CpGFrequency.pdf", plot = CpGFrequency, dpi = 300)

print("Plots Saved")

##### Final Model ####
#setting up predictors and response
predictors <- as.matrix(CpGmethylation %>% dplyr::select(-AgeY, -ID))
response <- as.matrix(CpGmethylation %>% dplyr::select(AgeY))

# fitting the model
glmnet_cv <- cv.glmnet(predictors, response, alpha = 0.05, type.measure = "mse")

#storing the lambda min 
full_lambda.min <- glmnet_cv$lambda.min

#running again for the lambda min
full_glmnet <- glmnet(predictors, response, alpha = 0.5, nlambda = 100)

# collecting up the non-zero CpGs and their coefficients
full_coef.lambda.min <- coef(full_glmnet, s = full_lambda.min) %>%
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "CpG") %>% 
  filter(s0 != 0)


#### Creating and saving useful outputs ####
write_csv(LOAOCV_output, paste0("OutputData/", runName ,"/LOAOCV_output.csv"))
write_csv(LOAOCV_predictions, paste0("OutputData/", runName, "/LOAOCV_predictions.csv"))
write_csv(variable_freq, paste0("OutputData/", runName, "/CpG_selection_freq.csv"))
write_csv(full_coef.lambda.min, paste0("OutputData/", runName, "/fullModelCoef.csv"))

print("Outputs written")