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
pairedData <- T

#UMI located in Read1 --> "R1"
#UMI located in Read1 and Read2 --> "R1 & R2"
UMIlocation <- "R1 & R2"

#length of the UMI
UMIlength <- 10

#length of th sequence
sequenceLength <- 251

#counts cutoff for initial UMI removal 
countsCutoff <- 0

#max UMI distance for UMI correction
UMIdistance <- 2

#max sequence distance for UMI correction
sequenceDistance <- 30

#filepaths
filepath1 <- "UMI in R1 and R2/BC_IGHFR1_S88_L001_R1_001.fastq.gz"
filepath2 <- "UMI in R1 and R2/BC_IGHFR1_S88_L001_R2_001.fastq.gz"

#outputs folder
outputsFolder <- "Outputs_Paired_R12_smallDS"

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
