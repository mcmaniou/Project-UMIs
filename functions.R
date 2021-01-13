groupingPairedR1 <-function(r1, counts, full, quality, full2, quality2, UMIlength){
  
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



groupingPairedR1R2 <-function(counts, full, quality, full2, quality2, r1){
  
  
  if ((counts == 1)){
    
    temp.result <- full[which(full$UMI12 == r1),]
    read1 <- temp.result$read
    read2 <- full2[which(full2$id == temp.result$id),read]
    
    quality.read1 <- as.character(quality[which(quality$id == temp.result$id), 1:(ncol(quality)-1)])
    quality.read1 <- as.numeric(quality.read1) + 33
    quality.read1 <- intToUtf8(quality.read1)
    
    quality.read2 <- as.character(quality2[which(quality2$id == temp.result$id), 1:(ncol(quality2)-1)])
    quality.read2 <- as.numeric(quality.read2) + 33
    quality.read2 <- intToUtf8(quality.read2)
    
    result <- data.table(UMI = r1,
                         read1 = read1, quality1 = quality.read1, 
                         read2 = read2, quality2 = quality.read2)
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
    
  }else{
    
    #reads with specific UMI
    grouping = full[which(full$UMI12 == r1), -"UMI12"]
    grouping2 = full2[which(full2$id %in% grouping$id), ]
    
    quality = quality[which(quality$id %in% grouping$id), ]
    quality2 = quality2[which(quality2$id %in% grouping2$id), ]
    
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
    
    result <- data.table(UMI = r1,
                         read1 = result1[1,1], quality1 = result1[1,2], 
                         read2 = result2[1,1], quality2 = result2[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  }
  return(result)
  
}



groupingSingle <-function(r1, counts, full, quality, UMIlength){
  
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



groupingFinalPairedR1 <-function(r1, full, quality, full2, quality2, first_consensus, UMIlength){
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



groupingFinalPairedR1R2 <- function(r1, full, quality, full2, quality2, first_consensus, UMIlength){
  if (nchar(r1) == 2*UMIlength){
    
    result <- first_consensus[which(first_consensus$UMI == r1),]
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
    
  }else{
    #reads with specific UMI
    grouping = full[str_detect(full$UMI12, r1), -"UMI12"]
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


one.run.calculationsFunction <- function(one.base){
  
  # setDT(exampleone.base)
  one.base = one.base[,.(mean = mean(V2), 
                       count = .N, 
                       perc_qual = mean(V2) / 93*100,
                       perc_count = .N/nrow(one.base)*100), by = list(V1),]
  
  one.base$criterion = rowMeans(one.base[,c("perc_qual", "perc_count")])
  
  #####!!! what happens if we have the same maximum criterion for different V1 !!!##########
 
  one.base <- one.base[which(one.base$criterion == max(one.base$criterion)), ]

  return(data.table(V1 = as.character(one.base[which.max(one.base$perc_qual),]$V1), 
                    mean = round(one.base[which.max(one.base$perc_qual),]$mean)))
  
  #return(data.table(V1 = as.character(one.base[which.max(one.base$criterion), ]$V1)[1], 
   #                 mean = round(one.base[which.max(one.base$criterion), ]$mean)[1]))
  
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



UMIcorrectionPairedR1 <- function(intermediate.table,first_consensus, sequenceDistance, UMIdistance){
  
  uniqueUMIs <- c()
  IDs_1 <- c()
  IDs_2 <- c()
  counts <- c()
  
  #while((nrow(intermediate.table)>1) & (intermediate.table[1,count] > 5)){
  while((nrow(intermediate.table)>1)){
    
    best <- first_consensus[which(first_consensus$UMI == intermediate.table$UMI[1]),]
    list.best <- best$UMI[1]
    list.counts <- intermediate.table$count[1]
    
    iterations <- c(2:nrow(intermediate.table))
    for (i in iterations){
      
      base_dist <- stringDist(c(best$UMI[1], intermediate.table$UMI[i]), method = "hamming") 
      
      if (base_dist <= UMIdistance){
        
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI[i]),read1]
        dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
        
        temp_read2 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI[i]),read2]
        dist2 <- stringDist(c(best$read2[1], temp_read2), method = "hamming") 
        
    
        if ((as.numeric(dist1) <= sequenceDistance) & (as.numeric(dist2) <= sequenceDistance)){
          list.best <- paste0(list.best,"|",intermediate.table$UMI[i])
          list.counts <- paste0(list.counts,"|",intermediate.table$count[i])
        }
      }
    }
    
    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2,intermediate.table$ID2[1])
    counts <- append(counts, list.counts)
    uniqueUMIs <- append(uniqueUMIs,list.best)
    intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI,as.character(list.best), negate = T), ]
    
  }
  
  if (!is.null(intermediate.table)){
    
    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2,intermediate.table$ID2[1])
    counts <- append(counts, intermediate.table$count[1])
    uniqueUMIs <- append(uniqueUMIs,intermediate.table$UMI[1])
    
  }
  
  newUMIs <- as.data.table(cbind(UMI = uniqueUMIs, ID1 = IDs_1, ID2 = IDs_2, Counts = counts))
  return(newUMIs)
}



UMIcorrectionPairedR1R2 <- function(intermediate.table,first_consensus, sequenceDistance, UMIdistance, UMIlength){
  
  uniqueUMIs <- c()
  IDs_1 <- c()
  IDs_2 <- c()
  counts <- c()
  
  #while((nrow(intermediate.table)>1) & (intermediate.table[1,count] > 5)){
  while((nrow(intermediate.table)>1)){
    
    best <- first_consensus[which(first_consensus$UMI == intermediate.table$UMI12[1]),]
    list.best.UMI <- best$UMI[1]
    list.counts <- intermediate.table$count[1]

    iterations <- c(2:nrow(intermediate.table))
    for (i in iterations){
      
      best.UMI1 <- substring(intermediate.table$UMI12[1],1,UMIlength)
      temp.UMI1 <- substring(intermediate.table$UMI12[i],1,UMIlength)
      best.UMI2 <- substring(intermediate.table$UMI12[1],UMIlength+1, 2*UMIlength)
      temp.UMI2 <- substring(intermediate.table$UMI12[i],UMIlength+1, 2*UMIlength)
      
      base_dist1 <- stringDist(c(best.UMI1, temp.UMI1), method = "hamming") 
      base_dist2 <- stringDist(c(best.UMI2, temp.UMI2), method = "hamming") 

      if ((base_dist1 <= UMIdistance) & (base_dist2 <= UMIdistance)){
       
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI12[i]),read1]
        dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
        
        temp_read2 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI12[i]),read2]
        dist2 <- stringDist(c(best$read2[1], temp_read2), method = "hamming") 
        
        if ((as.numeric(dist1) <= sequenceDistance) & (as.numeric(dist2) <= sequenceDistance)){
        
          list.best.UMI <- paste0(list.best.UMI,"|",intermediate.table$UMI12[i])
          list.counts <- paste0(list.counts,"|",intermediate.table$count[i])
          
        }
      }
    }
    
    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2,intermediate.table$ID2[1])
    counts <- append(counts,list.counts)
    uniqueUMIs <- append(uniqueUMIs,list.best.UMI)
  
    intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI12,as.character(list.best.UMI), negate = T), ]
    
  }
  
  if (!is.null(intermediate.table)){
    
    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    IDs_2 <- append(IDs_2,intermediate.table$ID2[1])
    counts <- append(counts, intermediate.table$count[1])
    uniqueUMIs <- append(uniqueUMIs,intermediate.table$UMI12[1])
    
  }
  
  newUMIs <- as.data.table(cbind(UMI = uniqueUMIs, ID1 = IDs_1, ID2 = IDs_2, Counts = counts))
  return(newUMIs)
}





UMIcorrectionSingle <- function(intermediate.table,first_consensus, sequenceDistance, UMIdistance){
  
  uniqueUMIs <- c()
  IDs_1 <- c()
  counts <- c()
  
  #while((nrow(intermediate.table)>1) & (intermediate.table[1,count] > 5)){
  while((nrow(intermediate.table)>1)){
    
    best <- first_consensus[which(first_consensus$UMI == intermediate.table$UMI[1]),]
    list.best <- best$UMI[1]
    list.counts <- intermediate.table$count[1]
    
    iterations <- c(2:nrow(intermediate.table))
    for (i in iterations){
      
      base_dist <- stringDist(c(best$UMI[1], intermediate.table$UMI[i]), method = "hamming") 
      
      if (base_dist <= UMIdistance){
        
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  intermediate.table$UMI[i]),read1]
        dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
            
        if (as.numeric(dist1) <= sequenceDistance) {
          list.best <- paste0(list.best,"|",intermediate.table$UMI[i])
          list.counts <- paste0(list.counts,"|",intermediate.table$count[i])
        }
      }
    }

    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    counts <- append(counts, list.counts)
    uniqueUMIs <- append(uniqueUMIs,list.best)
    intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI,as.character(list.best), negate = T), ]
    
  }
  
  if (!is.null(intermediate.table)){
    
    IDs_1 <- append(IDs_1,intermediate.table$ID1[1])
    counts <- append(counts, intermediate.table$count[1])
    uniqueUMIs <- append(uniqueUMIs,intermediate.table$UMI[1])
    
  }
  
  newUMIs <- as.data.table(cbind(UMI = uniqueUMIs, ID1 = IDs_1, Counts = counts))
  return(newUMIs)
}
