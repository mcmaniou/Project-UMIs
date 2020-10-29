groupingFunction <-function(r1, full, quality, full2, quality2){
  
  #print(r1)
  #reads with specific UMI
  # grouping = full[str_detect(full$UMI,as.character(r1)), ]
  
  grouping = full[which(full$UMI == r1), ]
  
  # grouping = as.data.frame(grouping)
  
  grouping2 <- full2[which(full2$id %in% grouping$id), ]
  
  quality = quality[which(quality$id %in% grouping$id), ]
  quality2 = quality2[which(quality2$id %in% grouping2$id), ]
  
  grouping = grouping[order(grouping$id), ]
  grouping2 = grouping2[order(grouping2$id), ]
  quality = quality[order(quality$id), ]
  quality2 = quality2[order(quality2$id), ]
  
  #File 1
  # grouping_q <- merge(grouping, quality, by ="id")
  
  grouping_q = cbind(grouping, quality[,-c("id")])
  
  #File 2
  grouping2$UMI <- grouping$UMI
  # grouping_q2 <- merge(grouping2, quality2, by="id")
  
  grouping_q2 = cbind(grouping2, quality2[,-c("id")])
  
  rm(grouping, grouping2, quality, quality2)
  
  result1 <- calculationsFunction(grouping_q) 
  result2 <- calculationsFunction(grouping_q2)
  
  result <- data.table(UMI = substr(r1,1,12),
                       read1 = result1[1,1], quality1 = result1[1,2], 
                       read2 = result2[1,1], quality2 = result2[1,2])
  
  colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  
  return(result)
}

groupingFinal <-function(r1, full, quality, full2, quality2, first_consensus){
  if (nchar(r1) == 12){
    
    result <- first_consensus[which(first_consensus$UMI == r1),]
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
    
  }else{
    #reads with specific UMI
    grouping = full[str_detect(full$UMI,as.character(r1)), ]
    
    # grouping = as.data.frame(grouping)
    
    grouping2 <- full2[which(full2$id %in% grouping$id), ]
    
    quality = quality[which(quality$id %in% grouping$id), ]
    quality2 = quality2[which(quality2$id %in% grouping2$id), ]
    
    grouping = grouping[order(grouping$id), ]
    grouping2 = grouping2[order(grouping2$id), ]
    quality = quality[order(quality$id), ]
    quality2 = quality2[order(quality2$id), ]
    
    #File 1
    # grouping_q <- merge(grouping, quality, by ="id")
    
    grouping_q = cbind(grouping, quality[,-c("id")])
    
    #File 2
    grouping2$UMI <- grouping$UMI
    # grouping_q2 <- merge(grouping2, quality2, by="id")
    
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

UMIcorrection <- function(test_full,first_consensus){
  
  uniqueUMIs <- c()

  while((nrow(test_full)>1) & (test_full[1,count] > 5)){

    best <- first_consensus[which(first_consensus$UMI == test_full$UMI[1]),]
    list.best <- best$UMI[1]
    
    iterations <- c(2:nrow(test_full))
    for (i in iterations){
      
      base_dist <- stringDist(c(best$UMI[1], test_full$UMI[i]), method = "hamming") 
      
      if (base_dist == 1){
        
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI[i]),read1]
        dist1 <- stringDist(c(best$read1[1], temp_read1), method = "hamming") 
        
        temp_read2 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI[i]),read2]
        dist2 <- stringDist(c(best$read2[1], temp_read2), method = "hamming") 
        
        #threshold <- 0.01 * (length(temp_read1) + length(temp_read2))
        #if (as.numeric(dist1 + dist2) < threshold){
        if ((as.numeric(dist1) < 4) & (as.numeric(dist2) < 4)){
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
