########## Loading libraries ########## 
library(tidyverse)
library(data.table)
library(ShortRead)
library(stringr)
library(Biostrings)
library(extrafont)
library(ggplot2)
library(hrbrthemes)

source("functions.R")

########## Inputs ##########
timeTable <- data.table(part = character(9),
                        time = character(9))

timeTable$part[1] <- "start time"
timeTable$time[1] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

fastqPath <- list.files(pattern = ".fastq", full = TRUE)
reads1 <- readFastq(fastqPath[1])
reads2 <- readFastq(fastqPath[2])

#table with the quality
qualityValues <- encoding(quality(reads1))

########## Data preparation ##########

#File 1
seq <- as.data.frame(sread(reads1))
ids <- as.data.frame(reads1@id)

full <- cbind(seq, ids)
names(full) <- c("seq", "id")

#the first 12 characters of the string are the UMI and the rest is the read
full$UMI <- substring(full$seq, 1, 12)
full$read <- substring(full$seq, 13, 251)
full <- separate(full,id, c("id1", "id2")," ", remove = T)
full <- select(full,read,id1,UMI)
colnames(full) <- c("read", "id", "UMI")

test_full<-full %>% 
  group_by(UMI) %>% 
  summarise(count = n())

quality <- as(quality(reads1), "matrix")
quality <- as.data.frame(quality)
quality$id <- full$id
quality<- as.matrix(quality)

#File 2
seq2 <- as.data.frame(sread(reads2))
ids2 <- as.data.frame(reads2@id)

full2 <- cbind(seq2,ids2)
names(full2) <- c("read", "id")
full2<- separate(full2,id, c("id1", "id2")," ",remove = T)
full2 <- select(full2,read, "id1")
colnames(full2) <- c("read", "id")

quality2 <- as(quality(reads2), "matrix")
quality2 <- as.data.frame(quality2)
quality2$id <- full2$id
quality2<- as.matrix(quality2)

########## Run approach 1 ##########
timeTable$part[2] <- "beforeApproach1"
timeTable$time[2] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

method <- "mean"
#consensus for each unique UMI
result_mean <- lapply(test_full$UMI, groupingFunction, full = full,
                      quality = quality, full2 = full2, quality2 = quality2, method = method)
result_mean <- bind_rows(result_mean)

result_sure <- result_mean

#UMI correction
timeTable$part[3] <- "beforeUMIcorrection1"
timeTable$time[3] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

newUMIs <- UMIcorrection(test_full,result_mean)

#consensus again for the new UMIs
timeTable$part[4] <- "beforeFinalConsensus1"
timeTable$time[4] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

consensus_mean <- data.table(UMI = character(),
                             read1 = character(), quality1 = character(), 
                             read2 = character(), quality2 = character())
iter <- c(1:length(newUMIs))
for (l in iter){
  input <- as.data.table(newUMIs[[l]])
  consensus_temp <- groupingNew(input, method)
  consensus_mean <- bind_rows(consensus_mean, consensus_temp)
  
}
#consensus_mean <- lapply(newUMIs, groupingNew, method = method)
#consensus_mean <- bind_rows(consensus_mean)

########## Run approach 2 ##########
timeTable$part[5] <- "beforeApproach2"
timeTable$time[5] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

method <- "multiply"

#consensus for each unique UMI
result_multiply <- lapply(test_full$UMI, groupingFunction, full = full,
                      quality = quality, full2 = full2, quality2 = quality2, method = method)
result_multiply <- bind_rows(result_multiply)

#UMI correction
timeTable$part[6] <- "beforeUMIcorrection2"
timeTable$time[6] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
 
newUMIsMul <- UMIcorrection(test_full,result_multiply)

#consensus again for the new UMIs
timeTable$part[7] <- "beforeFinalConsensus2"
timeTable$time[7] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

consensus_multiply <- data.table(UMI = character(),
                             read1 = character(), quality1 = character(), 
                             read2 = character(), quality2 = character())
iter <- c(1:length(newUMIsMul))
for (l in iter){
  input <- as.data.table(newUMIsMul[[l]])
  consensus_temp <- groupingNew(input, method)
  consensus_multiply <- bind_rows(consensus_multiply, consensus_temp)
  
}

#consensus_multiply <- lapply(newUMIsMul, groupingFunction,method = method)
#consensus_multiply <- bind_rows(consensus_multiply)


########## Outputs ##########
timeTable$part[8] <- "beforeOutputs"
timeTable$time[8] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

#approach 1-Mean
dir.create("Outputs Mean")

#File1
file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                   FastqQuality(consensus_mean$quality1),
                   BStringSet(consensus_mean$UMI))


fileSplit1 <- str_split(fastqPath[1],"\\.")
output <- paste0("Outputs Mean",fileSplit1[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

#File2
file <- ShortReadQ(DNAStringSet(consensus_mean$read2), 
                   FastqQuality(consensus_mean$quality2),
                   BStringSet(consensus_mean$UMI))

fileSplit2 <- str_split(fastqPath[2],"\\.")
output <- paste0("Outputs Mean",fileSplit2[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

#approach 2-Multiply
dir.create("Outputs Multiply")

#File1
file <- ShortReadQ(DNAStringSet(consensus_multiply$read1), 
                   FastqQuality(consensus_multiply$quality1),
                   BStringSet(consensus_multiply$UMI))

output <- paste0("Outputs Multiply",fileSplit1[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

#File2
file <- ShortReadQ(DNAStringSet(consensus_multiply$read2), 
                   FastqQuality(consensus_multiply$quality2),
                   BStringSet(consensus_multiply$UMI))

output <- paste0("Outputs Multiply",fileSplit2[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

remove(file, filepath, output, fileSplit1, fileSplit2)

timeTable$part[9] <- "end time"
timeTable$time[9] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

fwrite(timeTable, paste0(getwd(), "/time measurements.csv"),row.names = FALSE, sep = "\t", quote = FALSE)
