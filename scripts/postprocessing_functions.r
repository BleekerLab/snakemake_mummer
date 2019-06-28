# This file contains all functions required for postprocessing of the nucmer alignment data

# Reads the percentage aligned bases output files
get_aligned_perc_matrix <- function(filename = "results", Wsurpress = T){
  myfile <- filename
  if (!file.exists(myfile)){
    stop("your .tsv file was not found")
  }
  all_files <- list.files(".")
  tsv_files <- all_files[which(grepl(".tsv",all_files))] # Only tsv files
  if (length(tsv_files) > 1 & !Wsurpress){
    warning("multiple .tsv files, make sure you have the right one")
  }
  aligned_perc_matrix <- read.table(myfile,sep = "\t")
  #colnames(aligned_perc_matrix) <- 1:12
  return(aligned_perc_matrix)
}

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

# Locations of N in the DNA sequence
get_N_locations <- function(fastafile = opt$fasta, N_threshold = opt$n_threshold){
  my_seq <- readDNAStringSet(fastafile) 
  N_locations <- unlist(gregexpr("N",my_seq)) # Holy shit that are a lot of N, over 21% of the entire sequence
  if (N_locations[1] < 0){
    #print(paste0("A UNICORN! ",names(my_seq)," has NO N nucleotides!"))
    longasslist <- cbind(start = NA, end = NA, length = NA)
  } else {
    seq_mode <- 0
    longasslist <- c()
    for (i in 1:(length(N_locations)-1)){
      if(N_locations[i]+1==N_locations[i+1]){
        if(seq_mode!=1){
          mystart <- N_locations[i]
          seq_mode <- 1
        } else {
          #go on till you find the end
        }
      } else {
        if(seq_mode==1){
          myend <- N_locations[i]
          longasslist <- rbind(longasslist,c(start = mystart, end = myend)) # yay found a whole sequence
          seq_mode <- 2
        } else {
          longasslist <- rbind(longasslist,c(start = N_locations[i], end = N_locations[i])) # loner
          seq_mode <- 2
        }
      } 
    }
    if(seq_mode == 1){
      longasslist <- rbind(longasslist,c(start = mystart, end = N_locations[i+1]))
    }
    longasslist <- cbind(longasslist, length = longasslist[,2]-longasslist[,1]+1)
    #print(paste0(round(sum(longasslist[,3])/nchar(my_seq)*100,3),"% of ",names(my_seq)," is N"))
    #print(paste0(round(sum(longasslist[which(longasslist[,3]>=1000),3])/nchar(my_seq)*100,3),"% of ",names(my_seq)," is sequences of 1000 or more consecutive Ns")) #Still over 21%
    #write.table(longasslist, paste0(names(my_seq),"_N.tsv"))
    if (N_threshold>0){
      longasslist <- longasslist[-which(longasslist[,3]<N_threshold),]
    }
    if (!is.data.frame(longasslist)){
      longasslist <- cbind(start = longasslist[1], end = longasslist[2], length = longasslist[3])
    }
  }
  return(longasslist)
}

# NO LONGER NEEDED Generate table with on row = chr and col = SS and values are % aligned bases
# get_tab_align <- function(filenames = Out_align){
#   Tab_align <- c()
#   tmpcol <- c()
#   for (i in filenames){
#     tmp <- as.character(read.table(paste0("~/Wicher_comparitive_genomics/SS_alignments_nucmer/Output/",i))[,2])
#     Tab_align <- cbind(Tab_align, as.numeric(substr(tmp,1,nchar(tmp)-1)))
#     tmpcol <- c(tmpcol,substr(i,regexpr("_",i)+1,regexpr(".",i, fixed = T)-1))
#   }
#   colnames(Tab_align) <- tmpcol
#   rownames(Tab_align) <- paste0("H",1:12)
#   return(Tab_align)
# }
