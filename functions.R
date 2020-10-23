groupingFunction <-function(r1, full, quality, full2, quality2){
  
  #reads with specific UMI
  grouping = full[str_detect(full$UMI,as.character(r1)), ]
  grouping = as.data.frame(grouping)
  
  grouping2 <- full2[which(full2$id %in% grouping$id),]
  
  #File 1
  grouping_q<-merge(grouping, quality, by ="id")
  
  #File 2
  grouping2$UMI <- grouping$UMI
  grouping_q2 <- merge(grouping2, quality2, by="id")
  
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
    grouping = as.data.frame(grouping)
    
    grouping2 <- full2[which(full2$id %in% grouping$id),]
    
    #File 1
    grouping_q<-merge(grouping, quality, by ="id")
    
    #File 2
    grouping2$UMI <- grouping$UMI
    grouping_q2 <- merge(grouping2, quality2, by="id")
    
    result1 <- calculationsFunction(grouping_q) 
    result2 <- calculationsFunction(grouping_q2)
    
    result <- data.table(UMI = substr(r1,1,12),
                         read1 = result1[1,1], quality1 = result1[1,2], 
                         read2 = result2[1,1], quality2 = result2[1,2])
    
    colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  }
  
  return(result)
}

calculationsFunction <- function(grouping_q){
  
  cons <- matrix(nrow=length(grouping_q$read),ncol=nchar(grouping_q$read[1]))
  cons_corr <- matrix(nrow=2,ncol=nchar(grouping_q$read[1]))
  
  #loop runs for each sequence with the selected UMI
  #separates the sequence in letters
  loops <- c(1:length(grouping_q$read))
  for(i in loops){
    cons[i,]<- substring(grouping_q$read[i], 1:nchar(grouping_q$read[i]), 1:nchar(grouping_q$read[i]))
  }
  
  cons<- as.data.frame(cons)
  
  loops <- c(1: nchar(grouping_q$read[1]))
  
  #loop runs for each column/base
  for (y in loops){
    
    example <- cbind(as.character(cons[,y]), as.matrix(grouping_q[,(3+y)])) #(3+y) because the first two columns are "id", "UMI" and "read" 
    example <- as.data.frame(example)
    
    example$V2 <- as.numeric(example$V2)
  
    #groups per latter/base and finds the average quality
   
    table.q1 <- dplyr::count(example,V1) %>% 
      mutate(perc_count=n/length(cons[,1])*100)
    setDT(example)
    table.q2 <- example[,.(mean = mean(V2)),by=list(V1)]%>%  
       mutate(perc_qual=mean/93*100) 
    table.q <- merge(table.q1, table.q2 )%>% 
      mutate(criterion = rowMeans(select(.,c(perc_count, perc_qual)))) %>% 
      filter(criterion == max(criterion)) 
    
    #table.q <- example %>% group_by(V1) %>% summarise(count = n(), mean= mean(V2)) %>% 
     # mutate(perc_count=count/length(cons[,1])*100) %>%  #base/letter count percentage
    #  mutate(perc_qual=mean/93*100)   %>%  #93 quality values --> convert to percentage
     # mutate(criterion = rowMeans(select(.,c(perc_count, perc_qual)))) %>% #mean of the two percentages
    #  filter(criterion == max(criterion)) 
    
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


UMIcorrection <- function(test_full,first_consensus){
  
  uniqueUMIs <- data.table(UMI = character())
  #listNew <- list()
  #j <- 1
  
  while(nrow(test_full) > 0){
    
    maxCount <- max(test_full$count)
    if (maxCount < 6){
      break;
    }
    best <- test_full[which(test_full$count == maxCount ),]
    bestUMI <- as.character(best[1,1])
    bestCount <- as.numeric(best[1,2])
    bestread1 <- first_consensus[which(first_consensus$UMI == bestUMI), read1]
    bestread2 <- first_consensus[which(first_consensus$UMI == bestUMI), read2]
    list.best <- bestUMI
    #DT_UMI <- data.table( UMI = bestUMI, 
     #                     read1 = bestread1, 
      #                    quality1 = first_consensus[which(first_consensus$UMI == bestUMI), quality1],
       #                   read2 = bestread2,
        #                  quality2 = first_consensus[which(first_consensus$UMI == bestUMI), quality2],
         #                 counts = maxCount)
    
    iterations <- c(1:nrow(test_full))
    for (i in iterations){
      
      base_dist <- stringDist(c(as.character(bestUMI), test_full$UMI[i]), method = "hamming") 
      
      if (base_dist == 1){
        
        temp_read1 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI[i]),read1]
        dist1 <- stringDist(c(bestread1, temp_read1), method = "hamming") 
        
        temp_read2 <- first_consensus[which(first_consensus$UMI ==  test_full$UMI[i]),read2]
        dist2 <- stringDist(c(bestread2, temp_read2), method = "hamming") 
        
        threshold <- 0.05 * (length(temp_read1) + length(temp_read2))
        if (as.numeric(dist1 + dist2) < threshold){
         # DT_temp <- data.table( UMI = test_full$UMI[i], 
          #                       read1 = temp_read1, 
           #                      quality1 = first_consensus[which(first_consensus$UMI == test_full$UMI[i]), quality1],
            #                     read2 = temp_read2,
             #                    quality2 = first_consensus[which(first_consensus$UMI == test_full$UMI[i]), quality2],
              #                   counts = test_full$count[i])
          #DT_UMI <- bind_rows(DT_UMI, DT_temp)
          list.best <- paste0(list.best,"|",test_full$UMI[i])
        }
      }
    }
    
    list.new <-  data.table(UMI = character(1))
    list.new$UMI[1] <- as.character(list.best)
    uniqueUMIs <- bind_rows(uniqueUMIs, list.new)
    #listNew[[j]] <- as.data.table(DT_UMI)
    #j <- j+1
    test_full <- test_full[str_detect(test_full$UMI,as.character(list.best), negate = T), ]
    
  }
  
  return(uniqueUMIs)
}
