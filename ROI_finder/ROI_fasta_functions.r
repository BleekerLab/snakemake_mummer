# Functions for the ROI fasta script
library(OutlierDetection)
library(Biostrings)
library(seqinr)


# Get unique pieces of the chromosome names
get_chr_identifiers <- function(aligned_perc_matrix = aligned_perc_matrix){
  tmp <- colnames(aligned_perc_matrix)[1]
  for (h in 2:12){
    commonthing <- c()
    for (i in 1:length(unlist(strsplit(tmp,"")))){
      commonthing <- c(commonthing,unlist(strsplit(colnames(aligned_perc_matrix)[h],""))[i] == unlist(strsplit(tmp,""))[i])
    }
    tmp <- paste0(unlist(strsplit(colnames(aligned_perc_matrix)[h],""))[which(commonthing)], collapse = "")
  }
  chrID <- c()
  for (h in 1:12){
    uniquename <- c()
    for (i in 1:nchar(tmp)){
      uniquename <- c(uniquename, unlist(strsplit(colnames(aligned_perc_matrix)[h],""))[i] == unlist(strsplit(tmp,""))[i])
    }
    chrID <- c(chrID,paste0(unlist(strsplit(colnames(aligned_perc_matrix)[h],""))[-which(uniquename)], collapse = ""))
  }
  return(list(chrID = chrID, commonname = tmp))
}

# Assign super-scaffolds to chromosomes
assign_SStoChr <- function(accession = NA, plots = F){
  if (is.na(accession)){
    aligned_perc_matrix <- read.table("Input/results.tsv", sep = "\t")
  } else {
    aligned_perc_matrix <- read.table(paste0("Input/results",accession,".tsv"), sep = "\t") # get the aligned perc matrix
  }
  mylist <- get_chr_identifiers(aligned_perc_matrix = aligned_perc_matrix)
  chrID <- mylist$chrID
  commonname <- mylist$commonname
  colnames(aligned_perc_matrix) <- chrID
  aligned_perc_matrix <- aligned_perc_matrix[,order(as.numeric(colnames(aligned_perc_matrix)))] # fix its col names and order
  
  # Find the outliers in each of the entries
  Outliers <- c()
  for (i in 1:nrow(aligned_perc_matrix)){
    tmp <- OutlierDetection(t(aligned_perc_matrix[i,]), k = 0.05*nrow(t(aligned_perc_matrix)), rnames = T)
    if (!isEmpty(unname(tmp[2]))){
      Outliers <- c(Outliers,unname(tmp[[2]]))
      #print(paste0(i,": ",tmp[2]))
    } else {
      Outliers <- c(Outliers,0)
    }
  }
  
  # Find and plot the number of SS assigned to each chr
  freq_best_matches <- c()
  for (i in 0:12){
    freq_best_matches <- c(freq_best_matches, sum(Outliers==i))
  }
  names(freq_best_matches) <- 0:12
  if (plots){
    barplot(freq_best_matches, ylab = "matches", xlab = "chromosome", main = paste0(accession, " Super-Scaffolds aligned per chromsome"))
    for (i in 0:12){
      text(((i+1)*1.2)-0.5,2.5, labels = freq_best_matches[i+1])
    }
  }
  
  # return the matrix with the outliers in the next column for further processing
  aligned_perc_matrix <- cbind(aligned_perc_matrix,Outliers)
  return(aligned_perc_matrix)
}