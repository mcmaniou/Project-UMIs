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
start_time <- Sys.time()
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

uniqueUMIs <- UMIcorrection(full)

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
before1 <- Sys.time()
method <- "mean"
result_mean <- lapply(uniqueUMIs$UMI, groupingFunction, full = full,
                     quality = quality, full2 = full2, quality2 = quality2, method = method)
result_mean <- bind_rows(result_mean)

########## Run approach 2 ##########
before2 <- Sys.time()
method <- "multiply"
result_multiply <- lapply(uniqueUMIs$UMI, groupingFunction, full = full,
                          quality = quality, full2 = full2, quality2 = quality2, method = method)
result_multiply <- bind_rows(result_multiply)

########## Outputs ##########
beforeOutputs <- Sys.time()
#approach 1-Mean
dir.create("Outputs Mean")

#File1
file <- ShortReadQ(DNAStringSet(result_mean$read1), 
                    FastqQuality(result_mean$quality1),
                    BStringSet(result_mean$UMI))


fileSplit1 <- str_split(fastqPath[1],"\\.")
output <- paste0("Outputs Mean",fileSplit1[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

#File2
file <- ShortReadQ(DNAStringSet(result_mean$read2), 
                    FastqQuality(result_mean$quality2),
                    BStringSet(result_mean$UMI))

fileSplit2 <- str_split(fastqPath[2],"\\.")
output <- paste0("Outputs Mean",fileSplit2[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

#approach 2-Multiply
dir.create("Outputs Multiply")

#File1
file <- ShortReadQ(DNAStringSet(result_multiply$read1), 
                    FastqQuality(result_multiply$quality1),
                    BStringSet(result_multiply$UMI))

output <- paste0("Outputs Multiply",fileSplit1[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

#File2
file <- ShortReadQ(DNAStringSet(result_multiply$read2), 
                    FastqQuality(result_multiply$quality2),
                    BStringSet(result_multiply$UMI))

output <- paste0("Outputs Multiply",fileSplit2[[1]][2], "_corrected.fastq.gz")
file.create(output)
filepath <- paste0(output)
writeFastq(file, filepath, mode = "a")

remove(file, filepath, output, fileSplit1, fileSplit2)

########## Comparison ##########

#difference
base_difference <- data.table( diff1 = numeric(nrow(result_mean)),
                               diff2 = numeric(nrow(result_mean)))

iterations <- c(1:nrow(result_mean))
for (i in iterations){
  
  base_difference[i,1] <- stringDist(c(result_mean$read1[i],result_multiply$read1[i]), method = "hamming")
  base_difference[i,2] <- stringDist(c(result_mean$read2[i],result_multiply$read2[i]), method = "hamming")
  
}


fwrite(base_difference, "base difference.csv", row.names = FALSE, sep = "\t", quote = FALSE)

#font_import() 
loadfonts(device = "win")
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12, base_family = 'Times New Roman')+ 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

data.plot <- data.table(difference = c(base_difference$diff1, base_difference$diff2),
                         Type = c(rep_len("reads 1",nrow(base_difference)),rep_len("reads 2",nrow(base_difference))),
                         stringsAsFactors = F)

png(filename = paste(getwd(), "/Base difference.png", sep = ""), 
    width = 900, height = 700)

p <- ggplot(data=data.plot, aes(x=difference, group=Type, fill=Type)) +
  geom_density(adjust=1.5, alpha=.4) +
  ggtitle("Base difference") +
  xlab("Number of bases")+
  ylab("Density")+
  theme_ipsum()+
  theme(
    plot.title = element_text(size=16),
    axis.title.x = element_text(size = 15, vjust = 0.5,hjust = 0.5),
    axis.title.y = element_text(size = 15, vjust = 1.5,hjust = 0.5),
    legend.text = element_text(size=14),
    legend.title = element_text(size=15)
  )

print(p)
dev.off()

end_time <- Sys.time()
