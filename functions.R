groupingFunction <-function(r1, counts, full, quality, full2, quality2, UMIlength){
  
  #print(r1)
  
  if ((counts == 1)){
    
      temp.result <- full[which(full$UMI == r1),]
      read1 <- temp.result$read
      read2 <- full2[which(full2$id == temp.result$id),read]
      
      quality.read1 <- as.character(quality[which(quality$id == temp.result$id), 1:(ncol(quality)-1)])
      quality.read1 <- as.numeric(quality.read1) + 33
      quality.read1 <- intToUtf8(quality.read1)
      
      quality.read2 <- as.character(quality2[which(quality2$id == temp.result$id), 1:(ncol(quality2)-1)])
      quality.read2 <- as.numeric(quality.read2) + 33
      quality.read2 <- intToUtf8(quality.read2)
      
      result <- data.table(UMI = substr(r1,1,UMIlength),
                           read1 = read1, quality1 = quality.read1, 
                           read2 = read2, quality2 = quality.read2)
      
      colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  
  }else{
    
    #reads with specific UMI
    grouping = full[which(full$UMI == r1), ]
    grouping2 <- full2[which(full2$id %in% grouping$id), ]
    
    quality = quality[which(quality$id %in% grouping$id), ]
    quality2 = quality2[which(quality2$id %in% grouping2$id), ]
    
    grouping = grouping[order(grouping$id), ]
    grouping2 = grouping2[order(grouping2$id), ]
    quality = quality[order(quality$id), ]
    quality2 = quality2[order(quality2$id), ]
    
    #File 1
    grouping_q = cbind(grouping, quality[,-c("id")])
    
    #File 2
    grouping2$UMI <- grouping$UMI
    grouping_q2 = cbind(grouping2, quality2[,-c("id")])
    
    rm(grouping, grouping2, quality, quality2)
    
    result1 <- calculationsFunction(grouping_q) 
    result2 <- calculationsFunction(grouping_q2)
    
    result <- data.table(UMI = substr(r1,1,UMIlength),
                         read1 = result1[1,1], quality1 = result1[1,2], 
                         read2 = result2[1,1], quality2 = result2[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
    
  }
  
  return(result)
}



groupingPairedR1R2 <-function(ids, counts, full, quality, full2, quality2, UMI12){
  
  
  if ((counts == 1)){
    
    temp.result <- full[which(full$id == ids),]
    read1 <- temp.result$read
    read2 <- full2[which(full2$id == temp.result$id),read]
    
    quality.read1 <- as.character(quality[which(quality$id == temp.result$id), 1:(ncol(quality)-1)])
    quality.read1 <- as.numeric(quality.read1) + 33
    quality.read1 <- intToUtf8(quality.read1)
    
    quality.read2 <- as.character(quality2[which(quality2$id == temp.result$id), 1:(ncol(quality2)-1)])
    quality.read2 <- as.numeric(quality.read2) + 33
    quality.read2 <- intToUtf8(quality.read2)
    
    result <- data.table(UMI = UMI12,
                         read1 = read1, quality1 = quality.read1, 
                         read2 = read2, quality2 = quality.read2)
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
    
  }else{
    
    #reads with specific UMI
    grouping = full[str_detect(full$id,ids), ]
    grouping2 = full2[str_detect(full2$id,ids), ]
    
    quality = quality[str_detect(quality$id,ids), ]
    quality2 = quality2[str_detect(quality2$id,ids), ]
    
    grouping = grouping[order(grouping$id), ]
    grouping2 = grouping2[order(grouping2$id), ]
    quality = quality[order(quality$id), ]
    quality2 = quality2[order(quality2$id), ]
    
    #File 1
    grouping_q = cbind(grouping, quality[,-c("id")])
    
    #File 2
    grouping_q2 = cbind(grouping2, quality2[,-c("id")])
    
    rm(grouping, grouping2, quality, quality2)
    
    result1 <- calculationsFunction(grouping_q) 
    result2 <- calculationsFunction(grouping_q2)
    
    result <- data.table(UMI = UMI12,
                         read1 = result1[1,1], quality1 = result1[1,2], 
                         read2 = result2[1,1], quality2 = result2[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  }
  return(result)
  
}



groupingFunctionSingle <-function(r1, counts, full, quality, UMIlength){
  
  if ((counts == 1)){
    
    temp.result <- full[which(full$UMI == r1),]
    read1 <- temp.result$read
    
    quality.read1 <- as.character(quality[which(quality$id == temp.result$id), 1:(ncol(quality)-1)])
    quality.read1 <- as.numeric(quality.read1) + 33
    quality.read1 <- intToUtf8(quality.read1)
    
    result <- data.table(UMI = substr(r1,1,UMIlength),
                         read1 = read1, quality1 = quality.read1)
    
    colnames(result) <- c("UMI" , "read1", "quality1")
    
  }else{
  
    #print(r1)
    
    #reads with specific UMI
    grouping = full[which(full$UMI == r1), ]
    
    quality = quality[which(quality$id %in% grouping$id), ]
    
    grouping = grouping[order(grouping$id), ]
    
    quality = quality[order(quality$id), ]
    
    #File 1
    grouping_q = cbind(grouping, quality[,-c("id")])
    
    rm(grouping, quality)
    
    result1 <- calculationsFunction(grouping_q) 
    
    result <- data.table(UMI = substr(r1,1,UMIlength),
                         read1 = result1[1,1], quality1 = result1[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1")
  }
  
  return(result)
}



groupingFinal <-function(r1, full, quality, full2, quality2, first_consensus, UMIlength){
  if (nchar(r1) == UMIlength){
    
    result <- first_consensus[which(first_consensus$UMI == r1),]
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
    
  }else{
    #reads with specific UMI
    grouping = full[str_detect(full$UMI,as.character(r1)), ]
    grouping2 <- full2[which(full2$id %in% grouping$id), ]
    
    quality = quality[which(quality$id %in% grouping$id), ]
    quality2 = quality2[which(quality2$id %in% grouping2$id), ]
    
    grouping = grouping[order(grouping$id), ]
    grouping2 = grouping2[order(grouping2$id), ]
    quality = quality[order(quality$id), ]
    quality2 = quality2[order(quality2$id), ]
    
    #File 1
    grouping_q = cbind(grouping, quality[,-c("id")])
    
    #File 2
    grouping2$UMI <- grouping$UMI
    
    grouping_q2 = cbind(grouping2, quality2[,-c("id")])
    
    rm(grouping, grouping2, quality, quality2)
    
    result1 <- calculationsFunction(grouping_q) 
    result2 <- calculationsFunction(grouping_q2)
    
    result <- data.table(UMI = substr(r1,1,12),
                         read1 = result1[1,1], quality1 = result1[1,2], 
                         read2 = result2[1,1], quality2 = result2[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  }
  
  return(result)
}

groupingFinalPairedR1R2 <- function(r1, full, quality, full2, quality2, first_consensus, UMIlength, test_full){
  if (nchar(r1) == 2*UMIlength){
    
    result <- first_consensus[which(first_consensus$UMI == r1),]
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
    
  }else{
    #reads with specific UMI
    newIDs <- test_full[which(test_full$UMI12 %in% r1),"IDs"]
    grouping = full[which(full$id %in% newIDs$IDs), ]
    grouping2 <- full2[which(full2$id %in% newIDs$IDs), ]
    
    quality = quality[which(quality$id %in% newIDs$IDs), ]
    quality2 = quality2[which(quality2$id %in% newIDs$IDs), ]
    
    grouping = grouping[order(grouping$id), ]
    grouping2 = grouping2[order(grouping2$id), ]
    quality = quality[order(quality$id), ]
    quality2 = quality2[order(quality2$id), ]
    
    #File 1
    grouping_q = cbind(grouping, quality[,-c("id")])
    
    #File 2
    grouping_q2 = cbind(grouping2, quality2[,-c("id")])
    
    rm(grouping, grouping2, quality, quality2)
    
    result1 <- calculationsFunction(grouping_q) 
    result2 <- calculationsFunction(grouping_q2)
    
    result <- data.table(UMI = substr(r1,1,2*UMIlength),
                         read1 = result1[1,1], quality1 = result1[1,2], 
                         read2 = result2[1,1], quality2 = result2[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  }
  
  return(result)
}




groupingFinalSingle <-function(r1, full, quality, first_consensus, UMIlength){
  if (nchar(r1) == UMIlength){
    
    result <- first_consensus[which(first_consensus$UMI == r1),]
    colnames(result) <- c("UMI" , "read1", "quality1")
    
  }else{
    #reads with specific UMI
    grouping = full[str_detect(full$UMI,as.character(r1)), ]
    
    quality = quality[which(quality$id %in% grouping$id), ]
    
    grouping = grouping[order(grouping$id), ]
    
    quality = quality[order(quality$id), ]
    
    #File 1
    grouping_q = cbind(grouping, quality[,-c("id")])

    rm(grouping, quality)
    
    result1 <- calculationsFunction(grouping_q) 
    
    result <- data.table(UMI = substr(r1,1,UMIlength),
                         read1 = result1[1,1], quality1 = result1[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1")
  }
  
  return(result)
}


one.run.calculationsFunction <- function(example){
  
  # setDT(example)
  example = example[,.(mean = mean(V2), 
                       count = .N, 
                       perc_qual = mean(V2) / 93*100,
                       perc_count = .N/nrow(example)*100), by = list(V1),]
  
  example$criterion = rowMeans(example[,c("perc_qual", "perc_count")])
  
  #####!!! what happens if we have the same maximum criterion for different V1 !!!##########
 
  example <- example[which(example$criterion == max(example$criterion)), ]

  return(data.table(V1 = as.character(example[which.max(example$perc_qual),]$V1), 
                    mean = round(example[which.max(example$perc_qual),]$mean)))
  
  #return(data.table(V1 = as.character(example[which.max(example$criterion), ]$V1)[1], 
   #                 mean = round(example[which.max(example$criterion), ]$mean)[1]))
  
}

calculationsFunction <- function(grouping_q){
  
  # list with all nts
  cons = str_split(grouping_q$read, pattern = "", simplify = TRUE)
  cons = as.list(as.data.table(cons))
  
  # list with all qualities
  grouping_q = as.list(grouping_q[,4:ncol(grouping_q)])
  
  # merge lists 
  grouping_q = base::Map(data.table, cons, grouping_q)
  
  rm(cons)
  
  ##
  cons_corr = lapply(grouping_q, one.run.calculationsFunction)
  cons_corr = rbindlist(cons_corr)
  
  #join again in one final sequence
  consensus = str_c(cons_corr$V1, collapse = "")
  meanQuality <- as.numeric(cons_corr$mean) + 33 
  meanQuality <- intToUtf8(meanQuality)
  result <- data.table( seq = consensus, qual = meanQuality)
  return(result)
  
}



UMIcorrection <- function(test_full,first_consensus, sequenceDistance, UMIdistance){
  
  uniqueUMIs <- c()

  #while((nrow(test_full)>1) & (test_full[1,count] > 5)){
  while((nrow(test_full)>1)){
    
    best <- first_consensus[which(first_consensus$UMI == test_full$UMI[1]),]
    list.best <- best$UMI[1]
    
    iterations <- c(2:nrow(test_full))
    for (i in iterations){
      
      base_dist <- stringDist(c(best$UMI[1], test_full$UMI[i]), method = "hamming") 
      
      if (base_dist <= UMIdistance){
        
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI[i]),read1]
        dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
        
        temp_read2 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI[i]),read2]
        dist2 <- stringDist(c(best$read2[1], temp_read2), method = "hamming") 
        
    
        if ((as.numeric(dist1) <= sequenceDistance) & (as.numeric(dist2) <= sequenceDistance)){
          list.best <- paste0(list.best,"|",test_full$UMI[i])
        }
      }
    }
  
    uniqueUMIs <- append(uniqueUMIs,list.best)
    test_full <- test_full[str_detect(test_full$UMI,as.character(list.best), negate = T), ]
    
  }
  
  if (!is.null(test_full)){uniqueUMIs <- append(uniqueUMIs,test_full$UMI[1])}
  return(uniqueUMIs)
}



UMI12correction <- function(test_full,first_consensus, sequenceDistance, UMIdistance, UMIlength){
  
  uniqueUMIs <- c()
  
  #while((nrow(test_full)>1) & (test_full[1,count] > 5)){
  while((nrow(test_full)>1)){
    
    best <- first_consensus[which(first_consensus$UMI == test_full$UMI12[1]),]
    list.best.UMI <- best$UMI[1]

    iterations <- c(2:nrow(test_full))
    for (i in iterations){
      
      best.UMI1 <- substring(test_full$UMI12[1],1,UMIlength)
      temp.UMI1 <- substring(test_full$UMI12[i],1,UMIlength)
      best.UMI2 <- substring(test_full$UMI12[1],UMIlength+1, 2*UMIlength)
      temp.UMI2 <- substring(test_full$UMI12[i],UMIlength+1, 2*UMIlength)
      
      base_dist1 <- stringDist(c(best.UMI1, temp.UMI1), method = "hamming") 
      base_dist2 <- stringDist(c(best.UMI2, temp.UMI2), method = "hamming") 
      
      #print(base_dist1)
      #print(base_dist2)
      if ((base_dist1 <= UMIdistance) & (base_dist2 <= UMIdistance)){
       
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI12[i]),read1]
        dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
        
        temp_read2 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI12[i]),read2]
        dist2 <- stringDist(c(best$read2[1], temp_read2), method = "hamming") 
        print(list.best.UMI)
        print(i)
        print(dist1)
        print(dist2)
        
        if ((as.numeric(dist1) <= sequenceDistance) & (as.numeric(dist2) <= sequenceDistance)){
        
          list.best.UMI <- paste0(list.best.UMI,"|",test_full$UMI12[i])
        
        }
      }
    }
    
    uniqueUMIs <- append(uniqueUMIs,list.best.UMI)
  
    test_full <- test_full[str_detect(test_full$UMI12,as.character(list.best.UMI), negate = T), ]
    
  }
  
  if (!is.null(test_full)){uniqueUMIs <- append(uniqueUMIs,test_full$UMI12[1])}
  return(uniqueUMIs)
}





UMIcorrectionSingle <- function(test_full,first_consensus, sequenceDistance, UMIdistance){
  
  uniqueUMIs <- c()
  
  #while((nrow(test_full)>1) & (test_full[1,count] > 5)){
  while((nrow(test_full)>1)){
    
    best <- first_consensus[which(first_consensus$UMI == test_full$UMI[1]),]
    list.best <- best$UMI[1]
    
    iterations <- c(2:nrow(test_full))
    for (i in iterations){
      
      base_dist <- stringDist(c(best$UMI[1], test_full$UMI[i]), method = "hamming") 
      
      if (base_dist <= UMIdistance){
        
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI[i]),read1]
        dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
            
        if (as.numeric(dist1) <= sequenceDistance) {
          list.best <- paste0(list.best,"|",test_full$UMI[i])
        }
      }
    }
    
    uniqueUMIs <- append(uniqueUMIs,list.best)
    test_full <- test_full[str_detect(test_full$UMI,as.character(list.best), negate = T), ]
    
  }
  
  if (!is.null(test_full)){uniqueUMIs <- append(uniqueUMIs,test_full$UMI[1])}
  return(uniqueUMIs)
}
