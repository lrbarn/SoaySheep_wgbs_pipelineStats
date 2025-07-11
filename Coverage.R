#### Coverage Plots and Thresholds ####
####Packages####
library(tidyverse)

####Data####
# directory
dataDir <- "data/coverage"


genomeLength <- 2654047367

############################
# --- 1. Define your parameters ---
# These should be accurate to your actual file names and desired column names
coverage_thresholds <- c(4, 5, 6, 10, 15, 20, 25) # Add all the numeric thresholds you have files for
columns <- c("fileName", "NumBases") # Your desired column names
# If they are in a subfolder like 'data/coverage/', change this to 'data/coverage'

# --- 2. Create a list of all relevant files ---
# This glob pattern will match coverage_4x.txt, coverage_5x.txt, etc.
# The 'sprintf' makes the pattern flexible for single or double digit numbers
all_coverage_files <- list.files(
  path = dataDir,
  pattern = "^coverage_[0-9]+x\\.txt$", # Regex: starts with "coverage_", ends with ".txt", with numbers in between
  full.names = TRUE # Get the full path for read_csv
)

# --- 3. Read and process all files in one go using purrr::map_dfr() ---
# This is the core of the efficient solution
ImportCoverageData <- map_dfr(all_coverage_files, function(file_path) {
  
  # Extract the coverage threshold (e.g., "4x", "5x") from the filename
  # This regex captures the number followed by 'x'
  coverage_threshold <- str_extract(basename(file_path), "[0-9]+x")
  
  # Read the individual file
  df <- read_csv(
    file_path,
    col_names = FALSE, # Files have no headers
    show_col_types = FALSE, # Suppress messages about column types
    progress = FALSE # Suppress progress bar
  ) %>%
    # Assign column names
    set_names(columns) %>%
    # Add the coverage threshold column
    mutate(coverageThreshold = coverage_threshold)
  
  return(df)
})

### changing the coverage threshold to an ordinal factor
ImportCoverageData <- ImportCoverageData %>% 
  mutate(coverageThreshold = factor(coverageThreshold,
                                    levels = c("4x", "5x", "6x", "8x", "10x", "15x", "20x", "25x"),
                                    ordered = TRUE))

#### extract sample number and change to wide format ####
ImportCoverageData$sampleID <- as.factor(str_split_fixed(ImportCoverageData$fileName, "\\.", 3)[,1])
ImportCoverageData$perc_cov <- (ImportCoverageData$NumBases/genomeLength)*100

### wide data ###
coverageWide <- ImportCoverageData %>% 
  unique() %>% # there is some doubling up here, check the raw data
  select(-fileName) %>% 
  pivot_wider(names_from = coverageThreshold,
              values_from = NumBases,
              id_cols = sampleID)

# adding percentage 

coverageWide <- coverageWide %>% 
  mutate(perc_4x = (`4x` / genomeLength) * 100,
         perc_5x = (`5x` / genomeLength) * 100,
         perc_6x = (`6x` / genomeLength) * 100,
         perc_8x = (`8x` / genomeLength) * 100,
         perc_10x = (`10x` / genomeLength) * 100,
         perc_15x = (`15x` / genomeLength) * 100,
         perc_20x = (`20x` / genomeLength) * 100,
         perc_25x = (`25x` / genomeLength) * 100,)

#### summary stats for each threshold ####
summary_stats <- coverageWide %>% 
  summarise(
    # Use across() to apply functions to multiple columns
    # We select all columns EXCEPT SampleID, which is a character column
    across(.cols = -sampleID, # Selects all columns except 'SampleID'
           .fns = list(
             mean = ~ mean(.x, na.rm = TRUE),    # Calculate mean, ignoring NA values
             median = ~ median(.x, na.rm = TRUE), # Calculate median, ignoring NA values
             sd = ~ sd(.x, na.rm = TRUE),        # Calculate standard deviation, ignoring NA values
             min = ~ min(.x, na.rm = TRUE),      # Calculate minimum, ignoring NA values
             max = ~ max(.x, na.rm = TRUE)      # Calculate maximum, ignoring NA values              # Count non-NA values
           ),
           .names = "{.col}_{.fn}" # How to name the new columns (e.g., '4x_mean', '4x_sd')
    )
  )

sumstat <- ImportCoverageData %>% 
  group_by(coverageThreshold) %>% 
  summarise(mean = mean(NumBases),
            sd = sd(NumBases),
            min = min(NumBases),
            max = max(NumBases),
            mean_perc = mean(perc_cov),
            sd_perc = sd(perc_cov),
            min_perc = min(perc_cov),
            max_perc = max(perc_cov))


#### Plotting ####
ggplot(ImportCoverageData, aes(colour = coverageThreshold)) +
  geom_point(aes(x = coverageThreshold, y = perc_cov),alpha = 0.3) +
  geom_errorbar(data = sumstat, aes(x = coverageThreshold,
                                    ymin = mean_perc - sd_perc,
                                    ymax = mean_perc + sd_perc),
                width = 0.2) +
  geom_point(data = sumstat, aes(x = coverageThreshold, y = mean_perc), 
             size = 3, shape = 17, colour = "grey") +
  labs(x = "Coverage Threshold",
       y = "Genome Breadth (Percentage)") +
  theme_bw()

