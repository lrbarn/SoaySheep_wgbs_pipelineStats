# SoaySheep_wgbs_pipelineStats

this is the repository for looking into the quality of data from WGBS of Soay Sheep Samples

## Files
### Coverage.R
This collates the genome breadth for different coverage thresholds and creaates a plot

The data this script is expecting is a raw text file for each coverage threshold with the naming pattern coverage_4x.txt etc.
Within this there should be two columns of file name and raw bases covered.
The file name will follow the naming pattern sampleID.breadth_4x.txt etc
FOR EXAMPLE
the first line of coverage_20x.txt should be
001-1.breadth_20x.txt,200709239

### Pipeline_stats_LB.R
This collates some of the pipeline stats from multiQC data, samtools output, and BS conversion efficiency
Work in progress

The data expected here is a range
- fastQC shows how to extract the relevant fastQC data from Stanage
- Samtools Stats is expecting a text file called samtools_out.txt
    - this contains the selected data from the samtools_stats.txt files
    - how to extract this is in the Bash_cheatsheet
- BS Conversion Efficiency is expecting a table called BS_conversion_efficiency.txt
    - how to create this is in the Bash_Cheatsheet.txt
    - it is expecting three columns "max.conversion.failure", "min.conversion.efficiency", "sample.run.id"

### Bash_cheatsheet.txt
Useful bits of code for extracting the data from the Stanage HPC

### BASIS_SCRIPT.R
- this is a skeleton script for pulling together the raw .merged.bismakr.cov.gz files into methylKit
- then does some simple stats
- finally reformats the data into a large data frame to be used for elastic net regression in clock building

### MethylationAnalysis
- for running on stanage
- first 7 samples to make a LOOCV clock
- pdfs of graphs

## Directory Structure
- data is in a /data directory
  - coverage files are within in their own directory called /coverage
- scripts and bash_cheatsheet.txt etc all in the main working directory

# Miro Board Flowchart
https://miro.com/app/board/uXjVIgK9yZw=/

*Key* 
square = data
circle = script
rounded square = location
speech bubble = explanation
