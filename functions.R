groupingFunction <-function(r1, full, quality, full2, quality2, method){
  
  #reads with specific UMI
  grouping = full[str_detect(full$UMI,r1), ]
  grouping = as.data.frame(grouping)
  
  grouping2 <- full2[which(full2$id %in% grouping$id),]
  
  #File 1
  grouping_q<-merge(grouping, quality, by ="id")
  
  #File 2
  grouping2$UMI <- grouping$UMI
  grouping_q2 <- merge(grouping2, quality2, by="id")
  
  result1 <- calculationsFunction(grouping_q, method) 
  result2 <- calculationsFunction(grouping_q2, method)
  
  result <- data.table(UMI = substr(r1,1,12),
                       read1 = result1[1,1], quality1 = result1[1,2], 
                       read2 = result2[1,1], quality2 = result2[1,2])
  
  colnames(result) <- c("UMI" , "read1", "quality1", "read2", "quality2")
  
  return(result)
}


calculationsFunction <- function(grouping_q, method){
  
  cons <- matrix(nrow=length(grouping_q$read),ncol=nchar(grouping_q$read[1]))
  cons_corr <- matrix(nrow=2,ncol=nchar(grouping_q$read[1]))
  
  #loop runs for each sequence with the selected UMI
  #separates the sequence in letters
  loops <- c(1:length(grouping_q$read))
  for(i in loops){
    cons[i,]<- substring(grouping_q$read[i], 1:nchar(grouping_q$read[i]), 1:nchar(grouping_q$read[i]))
  }
  
  cons<- as.data.frame(cons)
  
  cons_corr <- switch(method, 
                      mean = methodMean(cons, cons_corr, grouping_q),
                      multiply = methodMultiply(cons, cons_corr, grouping_q)
                      )
  
  #join again in one final sequence
  consensus = str_c(cons_corr[1,], collapse = "")
  meanQuality <- as.numeric(cons_corr[2,]) +33 
  meanQuality <- intToUtf8(meanQuality)
  result <- data.table( seq = consensus, qual = meanQuality)
  return(result)
  
}


methodMultiply <- function(cons, cons_corr, grouping_q){
  
  loops <- c(1: nchar(grouping_q$read[1]))
  #loop runs for each column/base
  for (y in loops){
    
    example <- cbind(as.character(cons[,y]), as.matrix(grouping_q[,(3+y)])) #(3+y) because the first two columns are "id", "UMI" and "read" 
    example <-as.data.frame(example)
    
    example$V2 <- as.numeric(example$V2)
    
    #groups per latter/base and finds the average quality
    table.q <- example %>% group_by(V1) %>% summarise(count = n(), mean= mean(V2)) %>% 
      mutate(base_accuracy = (1-(10^(-mean/10)))) %>%
      mutate(criterion = base_accuracy * count) %>%
      filter(criterion == max(criterion))
    
    cons_corr[1,][y]<-as.character(table.q$V1)
    cons_corr[2,][y]<-round(table.q$mean)
  }
  
  return(cons_corr)
}


methodMean <- function(cons, cons_corr, grouping_q){
  
  loops <- c(1: nchar(grouping_q$read[1]))
  #loop runs for each column/base
  for (y in loops){
    
    example <- cbind(as.character(cons[,y]), as.matrix(grouping_q[,(3+y)])) #(3+y) because the first two columns are "id", "UMI" and "read" 
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
  
  return(cons_corr)
}

UMIcorrection <- function(full){
  
  test_full<-full %>% 
    group_by(UMI) %>% 
    summarise(count = n())
  
  uniqueUMIs <- data.table(UMI = character())
  
  while(nrow(test_full) > 0){
    
    maxCount <- max(test_full$count)
    if (maxCount < 6){
      break;
    }
    best <- test_full[which(test_full$count == maxCount ),]
    bestUMI <- best[1,1]
    bestCount <- best[1,2]
    list.best <- bestUMI
    
    iterations <- c(1:nrow(test_full))
    for (i in iterations){
      base_dist <- stringDist(c(as.character(bestUMI), test_full$UMI[i]), method = "hamming") 
      if ((base_dist == 1) && ( (bestCount/test_full$count) > 2 )){
        list.best <- paste0(list.best,"|",test_full$UMI[i]) 
      }
    }
    
    list.new <-  data.table(UMI = character(1))
    list.new$UMI[1] <- as.character(list.best)
    uniqueUMIs <- bind_rows(uniqueUMIs, list.new)
    test_full <- test_full[str_detect(test_full$UMI,as.character(list.best), negate = T), ]
    
  }
  
  return(uniqueUMIs)
}
