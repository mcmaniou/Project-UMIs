########## Loading libraries ########## 
library(tidyverse)
library(data.table)
library(ShortRead)
library(stringr)
library(Biostrings)

source("functions.R")

########## Inputs ##########
countsCutoff <- 5
distanceCutoff <- 4

timeTable <- data.table(part = character(6),
                        time = character(6))

timeTable$part[1] <- "start time"
timeTable$time[1] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

fastqPath <- list.files(pattern = ".fastq", full = TRUE)
reads1 <- readFastq(fastqPath[1])
reads2 <- readFastq(fastqPath[2])

#table with the quality
# qualityValues <- encoding(quality(reads1))

########## Data preparation ##########

#File 1
seq <- as.data.table(sread(reads1))
ids <- as.data.table(reads1@id)

full <- cbind(seq, ids)
names(full) <- c("seq", "id")

#the first 12 characters of the string are the UMI and the rest is the read
full$UMI <- substring(full$seq, 1, 12)
full$read <- substring(full$seq, 13, 251)
full <- separate(full, id, c("id1", "id2"), " ", remove = T)
full <- select(full, read, id1, UMI)
colnames(full) <- c("read", "id", "UMI")

test_full <- full[,.(count = .N),by=UMI,]
test_full <- test_full[which(test_full$count > countsCutoff),]

quality <- as(quality(reads1), "matrix")

quality = as.data.table(quality)

quality = quality[,13:251]

quality$id <- full$id

#File 2
seq2 <- as.data.table(sread(reads2))
ids2 <- as.data.table(reads2@id)

full2 <- cbind(seq2,ids2)
names(full2) <- c("read", "id")
full2<- separate(full2,id, c("id1", "id2")," ",remove = T)
full2 <- select(full2,read, "id1")
colnames(full2) <- c("read", "id")

quality2 <- as(quality(reads2), "matrix")

quality2 = as.data.table(quality2)

quality2$id <- full2$id

rm(ids, ids2, seq, seq2, reads1, reads2)

test_full = test_full[order(test_full$count, decreasing = TRUE), ]

########## First consensus ##########
timeTable$part[2] <- "beforeFirstConsensus"
timeTable$time[2] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

#test_full <- test_full[1:100,]
#consensus for each unique UMI
# result_mean <- lapply(test_full$UMI, groupingFunction, full = full,
#                       quality = quality, full2 = full2, quality2 = quality2)

result_mean = list()

for(i in test_full$UMI){
  
  result_mean[[i]] = groupingFunction(i, full, quality, full2, quality2)
  
}

result_mean <- bind_rows(result_mean)

########## UMI correction ##########
timeTable$part[3] <- "beforeUMIcorrection"
timeTable$time[3] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

newUMIs <- UMIcorrection(test_full,result_mean,distanceCutoff)

########## Final consensus ##########
timeTable$part[4] <- "beforeFinalConsensus"
timeTable$time[4] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

consensus_mean = list()

for(i in newUMIs){
  
  consensus_mean[[i]] = groupingFinal(i, full, quality, full2, quality2,result_mean)
  
}

consensus_mean <- bind_rows(consensus_mean)

#consensus_mean <- lapply(newUMIs$UMI, groupingFinal, full = full,
#                         quality = quality, full2 = full2, quality2 = quality2, first_consensus = result_mean)
#consensus_mean <- bind_rows(consensus_mean)

########## Outputs ##########
timeTable$part[5] <- "beforeOutputs"
timeTable$time[5] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

folder <- "Outputs"
dir.create(folder)

#File1
file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                   FastqQuality(consensus_mean$quality1),
                   BStringSet(consensus_mean$UMI))


fileSplit1 <- str_split(fastqPath[1],"\\.")
output <- paste0(folder,fileSplit1[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

#File2
file <- ShortReadQ(DNAStringSet(consensus_mean$read2), 
                   FastqQuality(consensus_mean$quality2),
                   BStringSet(consensus_mean$UMI))

fileSplit2 <- str_split(fastqPath[2],"\\.")
output <- paste0(folder,fileSplit2[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

remove(file, filepath, output, fileSplit1, fileSplit2)

timeTable$part[6] <- "end time"
timeTable$time[6] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")

fwrite(timeTable, paste0(getwd(), "/time measurements_10_3.csv"),row.names = FALSE, sep = "\t", quote = FALSE)
