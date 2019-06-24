# This file contains all the functions to get from a sorted coords file to a matrix with all percentages of aligned bases

# Get all relevant coords files filenames
get_coords_files <- function(keyphrase = "sortedSuper"){
  coords_files <- list.files("~/Wicher_comparitive_genomics/02SS_alignments_nucmer/Input")
  coords_files <- coords_files[which(grepl(keyphrase,coords_files))] # Only superscaffold results
  return(coords_files)
}

# Reads the sorted coords file, fixes its formatting and performs the filter step
get_filtered_data <- function(filename = "sortedSuper-Scaffold_1000002-H12.coords", length_threshold = opt$length_threshold, identity_threshold = opt$identity_threshold, mode = 0){
  if (!is.numeric(length_threshold)){stop("length_threshold needs to be a number")}
  if (!is.numeric(identity_threshold)){stop("identity_threshold needs to be a number")}
  if (identity_threshold>100 | identity_threshold<0){stop("identity_threshold needs to be a value between 0 and 100")}
  
  queryname <- substr(filename,regexpr("coords/",filename)+nchar("coords/"),regexpr("_vs_",filename)-1)
  refname <- substr(filename,regexpr("_vs_",filename)+nchar("_vs_"),regexpr(".sorted",filename)-1)
#  print(filename)
  # Read and fix
  if (length(scan(filename, what = character(), skip = 5, quiet = T))>0){
    data <- read.table(filename, sep = "", skip = 5)
    data <- data[,-which(data[1,]=="|")]
    colnames(data) <- tolower(c("Start_ref","End_ref","Start_query","End_query","Length_1","Length_2","%_Identity","Ref_name","Query_name"))
    
    # need to make start < end
    data[which(data$start_query > data$end_query),c(3,4)] <- data[which(data$start_query > data$end_query),c(4,3)]
    
    # filter
    if (all(data$length_2 < length_threshold | data$`%_identity` < identity_threshold)){
      warning("No alignments passed filtering, returning 0% aligned bases. Consider choosing different filtering parameters.")
      filtered_data <- cbind(rbind(rep(0,7)),ref_name = refname, query_name = queryname)
    } else if (any(data$length_2 < length_threshold | data$`%identity` < identity_threshold)){
      filtered_data <- data[-which(data$length_2 < length_threshold | data$`%_identity` < identity_threshold),]
    } else {
      filtered_data <- data
    }
    rownames(filtered_data) <- 1:nrow(filtered_data)
    filtered_data$ref_name <- refname
    filtered_data$query_name <- queryname
    if (mode == 0){
      write.table(filtered_data, paste0("Temp/filtered",filtered_data$ref_name[1],"-H",filtered_data$query_name[1],".tsv"), row.names = F, sep = "\t")
    } else if (mode == 1){
      return(filtered_data)
    }
  } else if (length(scan(filename, what = character(), skip = 5, quiet = T))==0){
    warning("file is empty, returning NA")
    filtered_data <- cbind(rbind(rep(NA,7)),ref_name = refname, query_name = queryname)
    if (mode == 0){
      write.table(filtered_data, paste0("Temp/filtered",filename,".tsv"), row.names = F, sep = "\t")
    } else if (mode == 1){
      return(filtered_data)
    }
  } else {
    stop(paste0(filename," can either not be read or somehow has a negative number of entries. Either way, something is really wrong"))
  }
}

# merge data
merge_data <- function(filename = "Temp/filteredSuper-Scaffold_1000002-H12.tsv",mode = 0, data = filtered_data, maxgap = 0){
  if (mode == 0){
    data <- read.table(filename, sep = "\t", header = T)
  }
  tbmerged <- data[,c(8,9,3,4)] # get just the relevant data

  
  if (all(is.na(tbmerged)[3:4])){
    merged <- tbmerged
  } else if (all((tbmerged == 0)[3:4])){
    merged <- tbmerged
  } else if (all(tbmerged[1,3:4] > 0)){
    tbmerged[,2] <- as.character(tbmerged[,2])
    tbmerged[,1] <- as.character(tbmerged[,1])
    while (T){
      merged <- c()
      i <- 1
      while (i<(nrow(tbmerged))){ # for loop which we can quickly skip through
        next_one <- which(tbmerged[(i+1):nrow(tbmerged),3]<=(tbmerged[i,4]+maxgap))
        if(length(next_one)>0){ # if true then that means we have overlap
          tmp <- which(tbmerged[,4]==max(c(tbmerged[next_one+i,4],0)))[1]
          merged <- rbind(merged,c(tbmerged[i,1],tbmerged[i,2],tbmerged[i,3],tbmerged[tmp,4]))
          i <- i+max(next_one)+1 # skip to the next not-overlapping sequence
        } else { # this sequence is not overlapping with others
          merged <- rbind(merged,c(tbmerged[i,1],tbmerged[i,2],tbmerged[i,3],tbmerged[i,4]))
          i <- i+1
        }
      }
      merged <- rbind(merged,c(tbmerged[nrow(tbmerged),1],tbmerged[nrow(tbmerged),2],tbmerged[nrow(tbmerged),3],tbmerged[nrow(tbmerged),4]))
      merged <- as.data.frame(merged, stringsAsFactors = F)
      merged[,3] <- as.numeric(merged[,3])
      merged[,4] <- as.numeric(merged[,4])
      if(nrow(tbmerged)==nrow(merged)){break}
      tbmerged <- merged
    }
    colnames(merged) <- c("ref_name","query_name","query_start","query_end")
  } else {stop("Data must be either NA or a numeric value >= 0")}
  if (mode == 0){
    write.table(merged, gsub(".tsv",".NR.tsv",filename), row.names = F, sep = "\t")
  } else {
    return(merged)
  }
}

# calculate the percentage of aligned bases
get_aligned_perc <- function(filename = "Temp/filteredSuper-Scaffold_1000002-H12.NR.tsv", fastafile = "Super-Scaffold_1000002.fasta", mode = 0, data = merged, out = opt$out, N_file = opt$nfile){
  tmp <- read.table(N_file, sep="\t", header = T)
  if (nrow(tmp)==0){
    totalN <- 0
  } else if (all(is.na(tmp)[1:3])){
    totalN <- 0
  } else {
    totalN <- sum(tmp[,3])
  }
  if (mode == 0){
    data <- read.table(filename, sep = "\t", header = T, stringsAsFactors = F)
  }
  if (all(is.na(data)[3:4])){
    aligned_perc <- cbind(ref_name = data[1],
                          query_name = data[2],
                          aligned_perc = NA)
  } else if (all((data == 0)[3:4])){
    aligned_perc <- cbind(ref_name = data[1],
                          query_name = data[2],
                          aligned_perc = 0)
  } else if (all(data[1,3:4] > 0)){
    data <- cbind(data,length = (data[,4]-data[,3])) # calculate length
    #  print(data)
    #  print(opt$filename)
    SS_seq <- readDNAStringSet(fastafile) 
    aligned_perc <- cbind(ref_name = data$ref_name[1],
                          query_name = data$query_name[1],
                          aligned_perc = sum(data$length)/(nchar(SS_seq)-totalN)*100) # length divided by total SS length *100 is % aligned bases
  } else {stop("Data must be either NA or a numeric value >= 0")}
  write.table(aligned_perc, out, sep = "\t", row.names = FALSE)
}

# filter the alignments that contain Ns
filter_N_entries <- function(data = filtered_data, N_file = opt$nfile){
  tmp <- read.table(N_file, sep="\t", header = T)
  if (nrow(tmp)==0){
    N_filtered_data <- data
  } else if (all(is.na(tmp)[1:3])){
    N_filtered_data <- data
  } else {
    N_containing_rows <- c()
    for (i in 1:nrow(tmp)){
      if (any(data[,3]<tmp[i,1] & data[,4]>tmp[i,2])){
        N_containing_rows <- c(N_containing_rows,which(data[,3]<tmp[i,1] & data[,4]>tmp[i,2]))
      }
    }
    if (is.null(N_containing_rows)){
      N_filtered_data <- data
    } else {
      N_filtered_data <- data[-N_containing_rows,]
    }
  }
  return(N_filtered_data)
}
