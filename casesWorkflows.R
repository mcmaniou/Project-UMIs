pairedR1 <- function(filepath1, filepath2, outputsFolder, 
                     UMIlength, UMIdistance, sequenceLength, sequenceDistance,
                     countsCutoff){
  
  #read input files
  reads1 <- readFastq(filepath1)
  reads2 <- readFastq(filepath2)
  
  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
  
  #separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength+1, sequenceLength)
  
  #File 2
  seq2 <- as.data.table(sread(reads2))
  ids2 <- as.data.table(reads2@id)
  
  full2 <- cbind(seq2,ids2)
  names(full2) <- c("read", "id")
   
  #keep intermediate information 
  partIDs <- as.data.table(cbind(UMI = full$UMI, ID1 = full$id, ID2 = full2$id))
  partIDs <- partIDs[!duplicated(UMI),]
  partIDs <- partIDs[order(partIDs$UMI, decreasing = TRUE), ]
  
  intermediate.table <- full[,.(count = .N),by=UMI,]
  intermediate.table <- intermediate.table[order(intermediate.table$UMI, decreasing = TRUE), ]
  intermediate.table$ID1 <- partIDs$ID1
  intermediate.table$ID2 <- partIDs$ID2
  intermediate.table <- intermediate.table[which(intermediate.table$count >= countsCutoff),]
  
  rm(partIDs)
  
  intermediate.table = intermediate.table[order(intermediate.table$count, decreasing = TRUE), ]
  
  #data preparation
  #file 1
  full <- separate(full, id, c("id1", "id2"), " ", remove = T)
  full <- select(full, read, id1, UMI)
  colnames(full) <- c("read", "id", "UMI")
  
  quality <- as(quality(reads1), "matrix")
  quality = as.data.table(quality)
  quality = quality[,(UMIlength+1):sequenceLength]
  quality$id <- full$id
  
  full2<- separate(full2,id, c("id1", "id2")," ",remove = T)
  full2 <- select(full2,read, "id1")
  colnames(full2) <- c("read", "id")
  
  quality2 <- as(quality(reads2), "matrix")
  quality2 = as.data.table(quality2)
  quality2$id <- full2$id
  
  rm(ids, ids2, seq, seq2, reads1, reads2)
  
  
  #first consensus
  result_mean = list()
  
  for(i in c(1:nrow(intermediate.table))){
    
    result_mean[[i]] = groupingPairedR1(intermediate.table$UMI[i], intermediate.table$count[i], full, quality, full2, quality2, UMIlength)
    
  }
  
  result_mean <- bind_rows(result_mean)
  
  #UMI correction
  newUMIs <- UMIcorrectionPairedR1(intermediate.table,result_mean,sequenceDistance, UMIdistance)
  rm(intermediate.table)
  #final consensus
  consensus_mean = list()
  
  for(i in newUMIs$UMI){
    
    consensus_mean[[i]] = groupingFinalPairedR1(i, full, quality, full2, quality2,result_mean, UMIlength)
    
  }
  
  consensus_mean <- bind_rows(consensus_mean)

  #produce Outputs 

  dir.create(outputsFolder)
  
  #File1
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(newUMIs$ID1," ",consensus_mean$UMI)))
  
  
  fileSplit <- as.data.table(str_split(filepath1,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  #File2
  file <- ShortReadQ(DNAStringSet(consensus_mean$read2), 
                     FastqQuality(consensus_mean$quality2),
                     BStringSet(paste0(newUMIs$ID2," ",consensus_mean$UMI)))
  
  fileSplit <- as.data.table(str_split(filepath2,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  
  output.csv <- as.data.table(cbind(UMI = consensus_mean$UMI,
                                    UMIs = newUMIs$UMI,
                                    counts = newUMIs$Counts,
                                    read1 = consensus_mean$read1,
                                    quality1 = consensus_mean$quality1,
                                    read2 = consensus_mean$read2,
                                    quality2 = consensus_mean$quality2,
                                    ID1 = newUMIs$ID1,
                                    ID2 = newUMIs$ID2))
  
  write.table(output.csv,paste0(outputsFolder,"/summary_table.csv"), sep ="\t", row.names = F)
  
  remove(file,output, fileSplit, output.csv)
  
}


pairedR1R2 <- function(filepath1, filepath2, outputsFolder, 
                       UMIlength, UMIdistance, sequenceLength, sequenceDistance,
                       countsCutoff){
  #read input files
  reads1 <- readFastq(filepath1)
  reads2 <- readFastq(filepath2)
  
  #reads1 <- reads1[1:1000]
  #reads2 <- reads2[1:1000]

  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
 
  #separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength+1, sequenceLength)
  
  #File 2
  seq2 <- as.data.table(sread(reads2))
  ids2 <- as.data.table(reads2@id)
  
  full2 <- cbind(seq2, ids2)
  names(full2) <- c("seq", "id")
  
  #separate UMI and read
  full2$UMI <- substring(full2$seq, 1, UMIlength)
  full2$read <- substring(full2$seq, UMIlength+1, sequenceLength)
  
  #keep the information of the read2 UMI, in the read1 also
  full$UMI12 <- paste0(full$UMI, full2$UMI)
  
  #keep intermediate information 
  partIDs <- as.data.table(cbind(UMI12 = full$UMI12, ID1 = full$id, ID2 = full2$id))
  
  intermediate.table <- partIDs[,.(count = .N), by = UMI12,]
  intermediate.table <- intermediate.table[which(intermediate.table$count >= countsCutoff),]
  
  partIDs <- partIDs[!duplicated(UMI12),]
  partIDs <- partIDs[order(partIDs$UMI12, decreasing = TRUE), ]
  
  intermediate.table <- intermediate.table[order(intermediate.table$UMI12, decreasing = TRUE), ]
  intermediate.table$ID1 <- partIDs$ID1
  intermediate.table$ID2 <- partIDs$ID2
  
  rm(partIDs)
  
  #data preparation
  full <- separate(full, id, c("id1", "id2"), " ", remove = T)
  full <- select(full, read, id1, UMI, UMI12)
  colnames(full) <- c("read", "id", "UMI", "UMI12")
    
  quality <- as(quality(reads1), "matrix")
  quality = as.data.table(quality)
  quality = quality[,(UMIlength+1):sequenceLength]
  quality$id <- full$id
  
  #File2
  full2 <- separate(full2, id, c("id1", "id2"), " ", remove = T)
  full2 <- select(full2, read, id1, UMI)
  colnames(full2) <- c("read", "id", "UMI")
  
  quality2 <- as(quality(reads2), "matrix")
  quality2 = as.data.table(quality2)
  quality2 = quality2[,(UMIlength+1):sequenceLength]
  quality2$id <- full2$id
  
  rm(ids, ids2, seq, seq2, reads1, reads2)
  
  #first consensus for each unique pair
  
  result_mean = list()
  
  for(i in c(1:nrow(intermediate.table))){
    
    result_mean[[i]] = groupingPairedR1R2(intermediate.table$count[i], full, quality, full2, quality2, intermediate.table$UMI12[i])
    
  }
  
  result_mean <- bind_rows(result_mean)
  
  #UMI correction
  #checks both UMI1 and UMI2
  newUMIs <- UMIcorrectionPairedR1R2(intermediate.table,result_mean,sequenceDistance, UMIdistance, UMIlength)
  rm(intermediate.table)
  
  #final consensus
  consensus_mean = list()
  
  for(i in c(1:nrow(newUMIs))){
    
    consensus_mean[[i]] = groupingFinalPairedR1R2(newUMIs$UMI[i], full, quality, full2, quality2,result_mean, UMIlength)
    
  }
  
  consensus_mean <- bind_rows(consensus_mean)
  
  #produce Outputs 
  
  dir.create(outputsFolder)
  
  #File1
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(newUMIs$ID1," ",substring(consensus_mean$UMI, 1,UMIlength))))
  
  
  fileSplit <- as.data.table(str_split(filepath1,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  #File2
  file <- ShortReadQ(DNAStringSet(consensus_mean$read2), 
                     FastqQuality(consensus_mean$quality2),
                     BStringSet(paste0(newUMIs$ID2," ",substring(consensus_mean$UMI, UMIlength+1, 2*UMIlength))))
  
  fileSplit <- as.data.table(str_split(filepath2,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  output.csv <- as.data.table(cbind(UMI = consensus_mean$UMI,
                                    UMIs = newUMIs$UMI,
                                    counts = newUMIs$Counts,
                                    read1 = consensus_mean$read1,
                                    quality1 = consensus_mean$quality1,
                                    read2 = consensus_mean$read2,
                                    quality2 = consensus_mean$quality2,
                                    ID1 = newUMIs$ID1,
                                    ID2 = newUMIs$ID2))
  
  write.table(output.csv,paste0(outputsFolder,"/summary_table.csv"), sep ="\t", row.names = F)
  
  remove(file, output, fileSplit, output.csv)
  
}

single <- function(filepath1, outputsFolder, 
                   UMIlength, UMIdistance, sequenceLength, sequenceDistance,
                   countsCutoff){
  
  #read input file
  reads1 <- readFastq(filepath1)
  
  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
  
  #separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength+1, sequenceLength)
  
  full <- select(full, read, id, UMI)
  colnames(full) <- c("read", "id", "UMI")
  
  
  #keep intermediate information 
  partIDs <- as.data.table(cbind(UMI = full$UMI, ID1 = full$id))
  partIDs <- partIDs[!duplicated(UMI),]
  partIDs <- partIDs[order(partIDs$UMI, decreasing = TRUE), ]
  
  intermediate.table <- full[,.(count = .N),by=UMI,]
  intermediate.table <- intermediate.table[order(intermediate.table$UMI, decreasing = TRUE), ]
  intermediate.table$ID1 <- partIDs$ID1
  intermediate.table <- intermediate.table[which(intermediate.table$count >= countsCutoff),]
  
  rm(partIDs)
  
  intermediate.table = intermediate.table[order(intermediate.table$count, decreasing = TRUE), ]
  
  #data preparation
  quality <- as(quality(reads1), "matrix")
  quality = as.data.table(quality)
  quality = quality[,(UMIlength+1):sequenceLength]
  quality$id <- full$id

  rm(ids, seq, reads1)
  
  #first consensus
  result_mean = list()
  
  for(i in c(1:nrow(intermediate.table))){
    
    result_mean[[i]] = groupingSingle(intermediate.table$UMI[i], intermediate.table$count[i], full, quality, UMIlength)
    
  }
  
  result_mean <- bind_rows(result_mean)
  
  #UMI correction
  newUMIs <- UMIcorrectionSingle(intermediate.table,result_mean,UMIdistance, sequenceDistance)
  rm(intermediate.table)
  
  #final consensus
  consensus_mean = list()
  
  for(i in newUMIs$UMI){
    
    consensus_mean[[i]] = groupingFinalSingle(i, full, quality,result_mean, UMIlength)
    
  }
  
  consensus_mean <- bind_rows(consensus_mean)
  
  #produce Outputs 
  
  dir.create(outputsFolder)
  
  #File1
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(newUMIs$ID1," ",consensus_mean$UMI)))
  
  
  fileSplit <- as.data.table(str_split(filepath1,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
 
  output.csv <- as.data.table(cbind(UMI = consensus_mean$UMI,
                                    UMIs = newUMIs$UMI,
                                    counts = newUMIs$Counts,
                                    read1 = consensus_mean$read1,
                                    quality1 = consensus_mean$quality1,
                                    ID1 = newUMIs$ID1))
  
  write.table(output.csv,paste0(outputsFolder,"/summary_table.csv"), sep ="\t", row.names = F)
  
  remove(file,output, fileSplit, output.csv)
  
}