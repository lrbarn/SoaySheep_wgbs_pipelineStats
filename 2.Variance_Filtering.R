#### 2.Variance_Filtering.R ####
## For filtering out sites that have the least variation in their methylation proportion
## follows on from 1.Methylation_Read_In.R

# ======= SPECIFY VARIANCE FILTER ===========
var_filter <- 0.01

print("Variance filter set to: ")
print(var_filter)

#### Packages ####
library(tidyverse)

#### Data ####
prop_methylated <- read_csv("OutputData/Prop_methylated.csv")

## Reformat to be able to create the summary stats
prop_methylatedLONG <- pivot_longer(prop_methylated, cols = 2:ncol(prop_methylated),
                                    names_to = "sampleID",
                                    values_to = "prop_meth")

#### Summary Stats ####
summary <- prop_methylatedLONG %>% 
  group_by(site) %>% 
  summarise(mean = mean(prop_meth),
            sd = sd(prop_meth),
            variance = var(prop_meth))

#### Plotting ####
variance_hist <- ggplot(summary, aes(x = variance)) +
  geom_histogram()

#### Filtering by Variance ####
highVarCpGs <- summary %>% 
  filter(variance >= var_filter)

print("Number of high var sites:")
print(nrow(highVarCpGs))

prop_methylated_filtered <- prop_methylated %>% 
  filter(site %in% highVarCpGs$site)

#### Saving outputs ####
write_csv(prop_methylated_filtered, "OutputData/prop_methylated_filtered.csv")
write_csv(summary, "OutputData/methylation_prop_summary.csv")
ggsave("varianceHistogram.pdf", plot = variance_hist, dpi = 300)
