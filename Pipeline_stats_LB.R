#### Multi QC analysis of Main data set ####

## Mostly working off Melanie's code but with some adjustments ###
## Data has been downloaded onto the local computer for now/ ease although this will likely need to be adjusted to be used straight on Stanage

#### Packages ####
## not sure what all of these do but they are at the top of melanie's script...
library(tidyverse)
library(stringr)
library(scico) # interfaace for colouring maps (pallets)
#library(rethinking)

#### Multi QC Data ####
SHORT_multiqc_fastqc_TRIMMED <- read_table("Documents/Sheffield/Projects/MLcode/data/SHORT_multiqc_fastqc_TRIMMED.txt") %>% 
  select(-6)

setwd("/Users/bip23lrb/Documents/Sheffield/Projects/MLcode")
## rename to match 
names(SHORT_multiqc_fastqc_TRIMMED) <- c("filename", "No_seq", "CG_percentage", "Dedups_perc", "Av_seq_length")

# extracting just the sample IDs
sample_ID_list <- str_split_fixed(SHORT_multiqc_fastqc_TRIMMED$filename, "_", 2)[,1]

# new column of the sample IDs
SHORT_multiqc_fastqc_TRIMMED$ID <- sample_ID_list

#at this point filtered to remove anything that has "R0" in the filename but it doesn't look like any in this dataset has this
# file is then saved at this point for use later

#write.table(SHORT_multiqc_fastqc_TRIMMED, "data/QC_trimmed.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#### Samtools Stats Data ####
## seeing what I can work with 
samtools_out <- read_delim("Documents/Sheffield/Projects/MLcode/data/samtools_out.txt", 
                           delim = ":", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)
names(samtools_out) <- c("variable", "value")
samtools_out$variable <- str_replace_all(samtools_out$variable, " ", "_")
samtools_out$value <- str_remove(samtools_out$value, ".samtools_stats.txt")

##THIS IS IN LONG FORMAT
# pivot table doesn't work so instead using vectors to rebuild the dataframe
sampleID <- samtools_out %>% 
  filter(variable == "sample") %>% 
  pull(value)

bases_mapped <- samtools_out %>% 
  filter(variable == "bases_mapped") %>% 
  pull(value)

error_rate <- samtools_out %>% 
  filter(variable == "error_rate") %>% 
  pull(value)

raw_total_sequences <- samtools_out %>% 
  filter(variable == "raw_total_sequences") %>% 
  pull(value)

reads_mapped_and_paired <- samtools_out %>% 
  filter(variable == "reads_mapped_and_paired") %>% 
  pull(value)

## Putting together in to a df
samtools_wide <- data.frame(sampleID = sampleID,
                            bases_mapped = bases_mapped,
                            error_rate = error_rate,
                            raw_total_sequences = raw_total_sequences,
                            reads_mapped_and_paired = reads_mapped_and_paired)

##THERE IS A MORE EFFICIENT WAY TO DO THIS

### plotting ###

#### BS Conversion Efficiency ####
BS_con <- read_table("Documents/Sheffield/Projects/SoaySheep_wgbs_pipelineStats/data/BS_conversion_efficiency.txt", col_names = FALSE)
names(BS_con) <- c("max.conversion.failure", "min.conversion.efficiency", "sample.run.id")

#splitting up sample name and lanes
BS_con <- BS_con %>% 
  mutate(sample = str_split_fixed(sample.run.id, "_", 2)[,1],
         lane = str_split_fixed(sample.run.id, "_", 3)[,3])

ggplot(BS_con, aes(x = min.conversion.efficiency)) +
  geom_bar() +
  theme_bw()

mean(BS_con$min.conversion.efficiency)
median(BS_con$min.conversion.efficiency)
