pairedR1 <- function(filepath1, filepath2, outputsFolder, 
                     UMIlength, UMIdistance, sequenceLength, sequenceDistance,
                     countsCutoff){
  
  #read input files
  reads1 <- readFastq(filepath1)
  reads2 <- readFastq(filepath2)
 
  #data preparation 
  
  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
  
  #separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength+1, sequenceLength)
  ID1 <- full$id[1]
  full <- separate(full, id, c("id1", "id2"), " ", remove = T)
  full <- select(full, read, id1, UMI)
  colnames(full) <- c("read", "id", "UMI")
  
  test_full <- full[,.(count = .N),by=UMI,]
  test_full <- test_full[which(test_full$count > countsCutoff),]
  
  quality <- as(quality(reads1), "matrix")
  
  quality = as.data.table(quality)
  
  quality = quality[,(UMIlength+1):sequenceLength]
  
  quality$id <- full$id
  
  #File 2
  seq2 <- as.data.table(sread(reads2))
  ids2 <- as.data.table(reads2@id)
  
  full2 <- cbind(seq2,ids2)
  names(full2) <- c("read", "id")
  ID2 <- full2$id[2]
  full2<- separate(full2,id, c("id1", "id2")," ",remove = T)
  full2 <- select(full2,read, "id1")
  colnames(full2) <- c("read", "id")
  
  quality2 <- as(quality(reads2), "matrix")
  
  quality2 = as.data.table(quality2)
  
  quality2$id <- full2$id
  
  rm(ids, ids2, seq, seq2, reads1, reads2)
  
  test_full = test_full[order(test_full$count, decreasing = TRUE), ]
  
  #first consensus
  result_mean = list()
  
  for(i in c(1:nrow(test_full))){
    
    result_mean[[i]] = groupingFunction(test_full$UMI[i], test_full$count[i], full, quality, full2, quality2, UMIlength)
    
  }
  
  result_mean <- bind_rows(result_mean)
  
  #UMI correction
  newUMIs <- UMIcorrection(test_full,result_mean,sequenceDistance, UMIdistance)
  
  #final consensus
  consensus_mean = list()
  
  for(i in newUMIs){
    
    consensus_mean[[i]] = groupingFinal(i, full, quality, full2, quality2,result_mean, UMIlength)
    
  }
  
  consensus_mean <- bind_rows(consensus_mean)

  #produce Outputs 

  dir.create(outputsFolder)
  
  #File1
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(ID1," ",consensus_mean$UMI)))
  
  
  fileSplit <- as.data.table(str_split(filepath1,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  #File2
  file <- ShortReadQ(DNAStringSet(consensus_mean$read2), 
                     FastqQuality(consensus_mean$quality2),
                     BStringSet(paste0(ID2," ",consensus_mean$UMI)))
  
  fileSplit <- as.data.table(str_split(filepath2,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  remove(file,output, fileSplit)
  
}


pairedR1R2 <- function(filepath1, filepath2, outputsFolder, 
                       UMIlength, UMIdistance, sequenceLength, sequenceDistance,
                       countsCutoff){
  #read input files
  reads1 <- readFastq(filepath1)
  reads2 <- readFastq(filepath2)
  
  #reads1 <- reads1[c(27601:28000,33501:33600,35501:36000)]
  #reads2 <- reads2[c(27601:28000,33501:33600,35501:36000)]
  #data preparation 
  
  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
 
  #separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength+1, sequenceLength)
  ID1 <- full$id[1]
  full <- separate(full, id, c("id1", "id2"), " ", remove = T)
  full <- select(full, read, id1, UMI)
  colnames(full) <- c("read", "id", "UMI")
    
  quality <- as(quality(reads1), "matrix")
  
  quality = as.data.table(quality)
  
  quality = quality[,(UMIlength+1):sequenceLength]
 
  quality$id <- full$id
  
  #File 2
  seq2 <- as.data.table(sread(reads2))
  ids2 <- as.data.table(reads2@id)
  
  full2 <- cbind(seq2, ids2)
  names(full2) <- c("seq", "id")

  #separate UMI and read
  full2$UMI <- substring(full2$seq, 1, UMIlength)
  full2$read <- substring(full2$seq, UMIlength+1, sequenceLength)
  ID2 <- full2$id[1]
  full2 <- separate(full2, id, c("id1", "id2"), " ", remove = T)
  full2 <- select(full2, read, id1, UMI)
  colnames(full2) <- c("read", "id", "UMI")
  
  quality2 <- as(quality(reads2), "matrix")
  
  quality2 = as.data.table(quality2)
  
  quality2 = quality2[,(UMIlength+1):sequenceLength]

  quality2$id <- full2$id
  
  rm(ids, ids2, seq, seq2, reads1, reads2)
  
  test_1 <- as.data.table(cbind(full$id, full$UMI))
  test_2 <- as.data.table(cbind(full2$id, full2$UMI))
  test_1 = test_1[order(test_1$V1, decreasing = TRUE), ]
  test_2 = test_2[order(test_2$V1, decreasing = TRUE), ]
  test_full <- as.data.table(cbind(test_1$V1, paste0(test_1$V2, test_2$V2)))
  colnames(test_full) <- c("ID","UMI12")
  
  uniquePairs <- test_full[,.(count = .N), by = UMI12,]
  uniquePairs <- uniquePairs[which(uniquePairs$count > countsCutoff),]
  
  test_full <- group_by(test_full,UMI12) %>%
    summarise(UMI12, IDs = paste(ID,collapse ="|"),)%>%
    unique()
  
  test_full <- inner_join(test_full,uniquePairs,by="UMI12")
  
  test_full = test_full[order(test_full$count, decreasing = TRUE), ]
  
  rm(uniquePairs)
  
  #first consensus for each unique pair
  
  result_mean = list()
  
  for(i in c(1:nrow(test_full))){
    
    result_mean[[i]] = groupingPairedR1R2(test_full$IDs[i], test_full$count[i], full, quality, full2, quality2, test_full$UMI12[i])
    
  }
  
  result_mean <- bind_rows(result_mean)
  
  #UMI correction
  #checks both UMI1 and UMI2
  test_full <- test_full[,1:2]
  newUMIs <- UMI12correction(test_full,result_mean,sequenceDistance, UMIdistance, UMIlength)
  
  #final consensus
  consensus_mean = list()
  
  for(i in newUMIs){
    
    consensus_mean[[i]] = groupingFinalPairedR1R2(i, full, quality, full2, quality2,result_mean, UMIlength, test_full)   
    
  }
  
  consensus_mean <- bind_rows(consensus_mean)
  
  #produce Outputs 
  
  dir.create(outputsFolder)
  
  #File1
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(ID1," ",consensus_mean$UMI)))
  
  
  fileSplit <- as.data.table(str_split(filepath1,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  #File2
  file <- ShortReadQ(DNAStringSet(consensus_mean$read2), 
                     FastqQuality(consensus_mean$quality2),
                     BStringSet(paste0(ID2," ",consensus_mean$UMI)))
  
  fileSplit <- as.data.table(str_split(filepath2,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
  
  remove(file,output, fileSplit)
  
}

single <- function(filepath1, outputsFolder, 
                   UMIlength, UMIdistance, sequenceLength, sequenceDistance,
                   countsCutoff){
  
  #read input file
  reads1 <- readFastq(filepath1)
  
  #data preparation 
  
  #File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)
  
  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")
  
  #separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength+1, sequenceLength)
  ID <- full$id[1]
  #full <- separate(full, id, c("id1", "id2"), " ", remove = T)
  full <- select(full, read, id, UMI)
  colnames(full) <- c("read", "id", "UMI")
  
  test_full <- full[,.(count = .N),by=UMI,]
  test_full <- test_full[which(test_full$count > countsCutoff),]
  
  quality <- as(quality(reads1), "matrix")
  
  quality = as.data.table(quality)
  
  quality = quality[,(UMIlength+1):sequenceLength]
  
  quality$id <- full$id
  
  
  rm(ids, seq, reads1)
  
  test_full = test_full[order(test_full$count, decreasing = TRUE), ]
  
  #first consensus
  result_mean = list()
  
  for(i in c(1:nrow(test_full))){
    
    result_mean[[i]] = groupingFunctionSingle(test_full$UMI[i], test_full$count[i], full, quality, UMIlength)
    
  }
  
  result_mean <- bind_rows(result_mean)
  
  #UMI correction
  newUMIs <- UMIcorrectionSingle(test_full,result_mean,UMIdistance, sequenceDistance)
  
  #final consensus
  consensus_mean = list()
  
  for(i in newUMIs){
    
    consensus_mean[[i]] = groupingFinalSingle(i, full, quality,result_mean, UMIlength)
    
  }
  
  consensus_mean <- bind_rows(consensus_mean)
  
  #produce Outputs 
  
  dir.create(outputsFolder)
  
  #File1
  file <- ShortReadQ(DNAStringSet(consensus_mean$read1), 
                     FastqQuality(consensus_mean$quality1),
                     BStringSet(paste0(ID," ",consensus_mean$UMI)))
  
  
  fileSplit <- as.data.table(str_split(filepath1,"\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)],"\\."))
  output <- paste0(outputsFolder,"/", fileSplit[1], "_corrected.fastq.gz")
  file.create(output)
  writeFastq(file, output, mode = "a")
 
  remove(file,output, fileSplit)
  
}