#### Site Variance ####
## working with the output of pilot7 which has the 7 random samples
## data was pulled together on stanage and the output just downloaded into the project area for testing
## goal here is to produce a script to look at and then filter by variance which will form a basis for a script to run on the full dataset on stanage

#### Packages ####
library(tidyverse)

#### Data ####
prop_methylated <- read_csv("data/Prop_methylated.csv")

#### Data Reformatting ####
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
## variance 
ggplot(summary, aes(x = variance)) +
  geom_histogram()

#### Filtering by variance ####
highVarCpGs_0.1 <- summary %>% 
  filter(variance >= 0.1)
highVarCpGs_0.1_vec <- c(highVarCpGs_0.1$site)


highVarCpGs_0.05 <- summary %>% 
  filter(variance >= 0.05)
highVarCpGs_0.05_vec <- c(highVarCpGs_0.05$site)

highVarCpGs_0.01 <- summary %>% 
  filter(variance >= 0.01)
highVarCpGs_0.01_vec <- c(highVarCpGs_0.01$site)

### Creating a new PropMethylated which has been filtered
prop_methylated_filtered <- prop_methylated %>% 
  filter(site %in% highVarCpGs_0.01$site)

write_csv(prop_methylated_fitlered, "OutputData/prop_methylated_filtered.csv")