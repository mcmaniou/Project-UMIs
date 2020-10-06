########## Loading libraries ########## 
library(tidyverse)
library(data.table)
library(ShortRead)
library(stringr)

########## Inputs ##########
fastqPath <- list.files(pattern = ".fastq", full = TRUE)
reads1 <- readFastq(fastqPath[1])
reads2 <- readFastq(fastqPath[2])
#fqFile <- FastqFile(fastqPath[1])

#table with the quality
qualityValues <- encoding(quality(reads1))

########## Data preparation ##########

#File 1
seq <- as.data.frame(sread(reads1))
qual <- quality(reads1)
qual <- as.data.frame(qual@quality)
ids <- as.data.frame(reads1@id)

full <- cbind(seq, qual, ids)
names(full) <- c("seq", "qual", "id")

#nchar(full$seq[1]) #251
#the first 12 characters of the string are the UMI and the rest is the read
full$UMI <- substring(full$seq, 1, 12)
full$read <- substring(full$seq, 13, 251)

test_full<-full %>% 
  group_by(UMI) %>% 
  summarise(count = n())

test_full<- test_full[which(test_full$count>5),]

quality <- as(quality(reads1), "matrix")
quality <- as.data.frame(quality)
quality$id <- ids

quality<- as.matrix(quality)
full <- separate(full,id, c("id1", "id2")," ", remove = F)

#File 2
seq2 <- as.data.frame(sread(reads2))
qual2 <- quality(reads2)
qual2 <- as.data.frame(qual2@quality)
ids2 <- as.data.frame(reads2@id)

full2 <- cbind(seq2, qual2, ids2)
names(full2) <- c("read", "qual", "id")

quality2 <- as(quality(reads2), "matrix")
quality2 <- as.data.frame(quality2)
quality2$id <- ids2

quality2<- as.matrix(quality2)
full2<- separate(full2,id, c("id1", "id2")," ",remove = F)

########## Functions ##########

calculations_mean <- function(grouping, grouping_q){
  
  cons <- matrix(nrow=length(grouping[,1]),ncol=nchar(grouping$read[1]))
  cons_corr <- matrix(nrow=2,ncol=nchar(grouping$read[1]))
  
  #loop runs for each sequence with the selected UMI
  #separates the sequence in letters
  loops <- c(1:length(grouping[,1]))
  for(i in loops){
    cons[i,]<- substring(grouping$read[i], 1:nchar(grouping$read[i]), 1:nchar(grouping$read[i]))
  }
  
  cons<- as.data.frame(cons)
  
  loops <- c(1: nchar(grouping$read[1]))
  #loop runs for each column/base
  for (y in loops){
    
    example <- cbind(as.character(cons[,y]), as.matrix(grouping_q[,(2+y)])) #(2+y) because the first two columns are "id" and "UMI" 
    example <-as.data.frame(example)
  
    example$V2 <- as.numeric(example$V2)

    #groups per latter/base and finds the average quality
    table.q <- example %>% group_by(V1) %>% summarise(count = n(), mean= mean(V2)) %>% 
      mutate(perc_count=count/length(cons[,1])*100) %>%  #base/letter count percentage
      mutate(perc_qual=mean/93*100)   %>%  #93 quality values --> convert to percentage
      mutate(criterion = rowMeans(select(.,c(perc_count, perc_qual)))) %>% #mean of the two percentages
      filter(criterion == max(criterion)) 
    
    cons_corr[1,][y]<-as.character(table.q$V1)
    cons_corr[2,][y]<-round(table.q$mean)
  }
  
  #join again in one final sequence
  consensus = str_c(cons_corr[1,], collapse = "")
  meanQuality <- as.numeric(cons_corr[2,]) +33 
  meanQuality <- intToUtf8(meanQuality)
  result <- data.table( seq = consensus, qual = meanQuality)
  return(result)
  
}

approach1_mean <-function(r1, full, quality, full2, quality2){
 
  #reads with specific UMI
  grouping = full[which(full$UMI == r1),]
  grouping = as.data.frame(grouping)
  
  grouping2 <- full2[which(full2$id1 %in% grouping$id1),]
  
  #File 1
  grouping_q<-merge(grouping[,c(3,6)], quality, by ="id")
  result1 <- calculations_mean(grouping, grouping_q)  
  
  #File 2
  grouping2$UMI <- grouping$UMI
  grouping_q2<-merge(grouping2[,c(3,6)], quality2, by="id")
  result2 <- calculations_mean(grouping2, grouping_q2)
  
  
  result <- data.table(UMI = r1, read1 = result1[1,1], quality1 = result1[1,2], read2 = result2[1,1], quality2 = result2[1,2])
  colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  
  return(result)
}

########## Run approach 1 ##########

result <- lapply(test_full$UMI, approach1_mean, full = full, quality = quality, full2 = full2, quality2 = quality2)

result <- bind_rows(result)

########## Outputs ##########
dir.create("Outputs")

#File1
file1 <- ShortReadQ(DNAStringSet(result$read1), 
                    FastqQuality(result$quality1),
                    BStringSet(result$UMI))

fileSplit1 <- str_split(fastqPath[1],"\\.")
output1 <- paste0("Outputs",fileSplit1[[1]][2], "_corrected.fastq.gz")
file.create(output1)
filepath1 <- paste0(output1)
writeFastq(file1, filepath1, mode = "a")

#File2
file2 <- ShortReadQ(DNAStringSet(result$read2), 
                    FastqQuality(result$quality2),
                    BStringSet(result$UMI))

fileSplit2 <- str_split(fastqPath[2],"\\.")
output2 <- paste0("Outputs",fileSplit2[[1]][2], "_corrected.fastq.gz")
file.create(output2)
filepath2 <- paste0(output2)
writeFastq(file2, filepath2, mode = "a")

