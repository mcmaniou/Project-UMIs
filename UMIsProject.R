########## Loading libraries ########## 
library(tidyverse)
library(data.table)
library(ShortRead)
library(stringr)
library(Biostrings)

source("casesWorkflows.R")
source("functions.R")

########## Inputs ##########

#type of data - paired or single
pairedData <- F

#UMI located in Read1 --> "R1"
#UMI located in Read1 and Read2 --> "R1 & R2"
UMIlocation <- "R1"

#length of the UMI
UMIlength <- 12

#length of th sequence
sequenceLength <- 251

#min read counts per UMI, for initial data cleaning
countsCutoff <- 5

#max UMI distance for UMI merging
UMIdistance <- 1

#max sequence distance for UMI correction
sequenceDistance <- 3

#filepaths
filepath1 <- "UMI in R1/BC01_1680_IGHFR1_S89_L001_R1_001.fastq.gz"
filepath2 <- "UMI in R1/BC01_1680_IGHFR1_S89_L001_R2_001.fastq.gz"

#outputs folder
outputsFolder <- "Outputs_Single"

########## Run the appropriate scenario ##########


if (pairedData & UMIlocation == "R1"){   #case 1 -- paired data and UMI only in Read1
  
  pairedR1(filepath1, filepath2, outputsFolder, 
           UMIlength, UMIdistance, sequenceLength, sequenceDistance, countsCutoff)
  
}else if (pairedData & UMIlocation == "R1 & R2"){   #case 2 -- paired data and UMI in Read1 and Read2
  
  pairedR1R2(filepath1, filepath2, outputsFolder, 
             UMIlength, UMIdistance, sequenceLength, sequenceDistance, countsCutoff)
  
}else if (!pairedData){  #case 3 -- single data

  single(filepath1, outputsFolder, 
         UMIlength, UMIdistance, sequenceLength, sequenceDistance, countsCutoff)
  
}  
