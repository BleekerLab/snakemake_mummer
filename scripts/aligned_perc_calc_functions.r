# This file contains all the functions to get from a sorted coords file to a matrix with all percentages of aligned bases
# The functions in this file are used by filter_data.r and ref-query_matching.r

# Get all relevant coords files filenames
# No longer used!
get_coords_files <- function(keyphrase = "sortedSuper"){
  coords_files <- list.files("~/Wicher_comparitive_genomics/02SS_alignments_nucmer/Input")
  coords_files <- coords_files[which(grepl(keyphrase,coords_files))] # Only superscaffold results
  return(coords_files)
}

# Reads the sorted coords file, fixes its formatting and performs the filter step
get_filtered_data <- function(filename = "sortedSuper-Scaffold_1000002-H12.coords", length_threshold = opt$length_threshold, identity_threshold = opt$identity_threshold, mode = 1){
  if (!is.numeric(length_threshold)){stop("length_threshold needs to be a number")}
  if (!is.numeric(identity_threshold)){stop("identity_threshold needs to be a number")}
  if (identity_threshold>100 | identity_threshold<0){stop("identity_threshold needs to be a value between 0 and 100")}
  
  # Get the names of the query and reference for future use
  queryname <- substr(filename,regexpr("coords/",filename)+nchar("coords/"),regexpr("_vs_",filename)-1)
  refname <- substr(filename,regexpr("_vs_",filename)+nchar("_vs_"),regexpr(".sorted",filename)-1)

  # Read and fix the table layout 
  # The .coords files contain 5 leading rows before the actual table starts, so we use skip = 5 to skip these and avoid formatting issues.
  # If, after skipping 5 rows we're left with an empty file then we have no alignments, most likely caused by too strict alignment parameters in nucmer.
  if (length(scan(filename, what = character(), skip = 5, quiet = T))>0){
    data <- read.table(filename, sep = "", skip = 5)
    data <- data[,-which(data[1,]=="|")]
    colnames(data) <- tolower(c("Start_ref","End_ref","Start_query","End_query","Length_1","Length_2","%_Identity","Ref_name","Query_name"))
    
    # need to make start < end (force forward orientation)
    data[which(data$start_query > data$end_query),c(3,4)] <- data[which(data$start_query > data$end_query),c(4,3)]
    data <- data[order(data[,3]),] # need to reorder the file
    rownames(data) <- 1:nrow(data)
    
    # filter based on length and identity score
    if (all(data$length_2 < length_threshold | data$`%_identity` < identity_threshold)){
      warning("No alignments passed filtering, returning 0% aligned bases. Consider choosing different filtering parameters.")
      filtered_data <- as.data.frame(cbind(rbind(rep(0,7)),ref_name = refname, query_name = queryname), stringsAsFactors = FALSE)
    } else if (any(data$length_2 < length_threshold | data$`%_identity` < identity_threshold)){
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
    filtered_data <- as.data.frame(cbind(rbind(rep(NA,7)),ref_name = refname, query_name = queryname), stringsAsFactors = FALSE)
    if (mode == 0){
      write.table(filtered_data, paste0("Temp/filtered",filename,".tsv"), row.names = F, sep = "\t")
    } else if (mode == 1){
      return(filtered_data)
    }
  } else {
    stop(paste0(filename," can either not be read or somehow has a negative number of entries. Either way, something is really wrong"))
  }
}

# filter the alignments that contain Ns
filter_N_entries <- function(data = filtered_data, N_file = opt$nfile, N_threshold_perc = opt$n_threshold_perc, mode = 1, outname = opt$out){
  tmp <- read.table(N_file, sep="\t", header = T)
  if (nrow(tmp)==0){ #empty tmp file means there are no N nucleotides in the sequence
    N_filtered_data <- data
  } else if (all(is.na(tmp)[1:3])){ #don't know how this can happen but it apparently can so yeah... needed to have a fix for that too
    N_filtered_data <- data
  } else { #not empty, not messed up, lets roll
    N_containing_rows <- c()
    for (i in 1:nrow(data)){
      N_amount <- 0
      # check if there is any overlap between the alignment and the sequences of N nucleotides in the sequence(N-sequence).
      # This looks complex but really it ain't.
      if (any(tmp[,1]>=data[i,3] & tmp[,1]<=data[i,4]) #does the beginning of any of the N-sequence fall within the alignment
           | any(tmp[,2]>=data[i,3] & tmp[,2]<=data[i,4])){ #does the end of any of the N-sequence fall within the alignment
        for (j in which((tmp[,1]>=data[i,3] & tmp[,1]<=data[i,4]) | (tmp[,2]>=data[i,3] & tmp[,2]<=data[i,4]))){
          if (tmp[j,1]>=data[i,3] & tmp[j,2]<=data[i,4]){ #the entirety of the N-sequence falls within the alignment
            N_amount <- N_amount + tmp[j,3] #add the length of the N-sequence to the N_amount for the current alignment
          } else if (tmp[j,1]<data[i,3] & tmp[j,2]<=data[i,4]){ #the beginning of the N-sequence is outside the alignment but the end is within
            N_amount <- N_amount + tmp[j,2]-data[i,3] #only add the piece that falls within the alignment to N_amount
          } else if (tmp[j,1]>=data[i,3] & tmp[j,2]>data[i,4]){ #the beginning of the N-sequence is within the alignment but the end is outside
            N_amount <- N_amount + data[i,4]-tmp[j,1] #only add the piece that falls within the alignment to N_amount
          }
        }
        
        # Check if the % of N nucleotides exceeds the threshold value
        if (100*N_amount/(data[i,4]-data[i,3]) > N_threshold_perc){
          N_containing_rows <- c(N_containing_rows,i)
        }
      } else if (any(tmp[,1]<data[i,3] & tmp[,2]>data[i,4])){
        N_containing_rows <- c(N_containing_rows,i)
      }
    }
    
    # Filter out the alignments that contained too much N
    if (is.null(N_containing_rows)){
      N_filtered_data <- data
    } else {
      N_filtered_data <- data[-N_containing_rows,]
    }
  }
  
  if (nrow(N_filtered_data)==0){
    warning("No alignments passed N_filtering, returning 0% aligned bases. Consider choosing different filtering parameters.")
    N_filtered_data <- as.data.frame(cbind(rbind(rep(0,7)), ref_name = data$ref_name[1], query_name = data$query_name[1]), stringsAsFactors = FALSE)
  }
  if (mode == 0){
    write.table(N_filtered_data, outname, row.names = F, sep = "\t")
  } else {
    return(N_filtered_data)
  }
}

# merge data
# Should be replaced by bedtools
merge_data <- function(filename = "Temp/filteredSuper-Scaffold_1000002-H12.tsv",mode = 1, data = filtered_data, maxgap = 0, datatype = 1, maxgap2 = 0){
  if (mode == 0){
    data <- read.table(filename, sep = "\t", header = T)
  }
  if (datatype == 1){
    tbmerged <- data[,c(8,9,3,4)] # get just the relevant query data
  } else if (datatype == 2){
    tbmerged <- data[,c(8,9,1,2)] # get just the relevant reference
  } else if (datatype == 3){
    tbmerged <- data[,c(8,9,3,4,1,2,6,5,7)] # get everything but still use query formatting
  } else if (datatype == 4){
    tbmerged <- data[,c(8,9,1:7)] # get everything but still use ref formatting
  } else if (datatype == 5){
    tbmerged <- data[,-5]
  } else {
    stop("datatype needs to be either 1,2,3,4, or 5")
  }
  

  #print("merge step 1")
  if (all(is.na(tbmerged)[3:4])){
    merged <- tbmerged
  } else if (all((tbmerged == 0)[3:4])){
    merged <- tbmerged
  } else if (nrow(tbmerged)==1){
    merged <- tbmerged
  } else if (all(tbmerged[1,3:4] > 0)){
    tbmerged[,2] <- as.character(tbmerged[,2])
    tbmerged[,1] <- as.character(tbmerged[,1])
    while (T){
      merged <- data.frame()
      i <- 1
      while (i<(nrow(tbmerged))){ # for loop which we can quickly skip through
        next_one <- which(tbmerged[(i+1):nrow(tbmerged),3]<=(tbmerged[i,4]+maxgap))
        if(length(next_one)>0){ # if true then that means we have overlap
          tmp <- which(tbmerged[,4]==max(c(tbmerged[next_one+i,4],tbmerged[i,4])))[1]
          merged <- rbind(merged,cbind(tbmerged[i,1],tbmerged[i,2],tbmerged[i,3],tbmerged[tmp,4]), stringsAsFactors = FALSE)
          i <- i+max(next_one)+1 # skip to the next not-overlapping sequence
          if (i==(nrow(tbmerged))){
            lastmerged <- FALSE
          } else {
            lastmerged <- TRUE
          }
        } else { # this sequence is not overlapping with others
          merged <- rbind(merged,cbind(tbmerged[i,1],tbmerged[i,2],tbmerged[i,3],tbmerged[i,4]), stringsAsFactors = FALSE)
          i <- i+1
          lastmerged <- FALSE
        }
      }
      if (!lastmerged){
        merged <- rbind(merged,cbind(tbmerged[nrow(tbmerged),1],tbmerged[nrow(tbmerged),2],tbmerged[nrow(tbmerged),3],tbmerged[nrow(tbmerged),4]), stringsAsFactors = FALSE)
      }
      merged <- as.data.frame(merged, stringsAsFactors = F)
      merged[,3] <- as.numeric(merged[,3])
      merged[,4] <- as.numeric(merged[,4])
      if(nrow(tbmerged)==nrow(merged) | nrow(merged) == 1){break}
      tbmerged <- merged
    }
    
    #print("merge step 2")
    if (maxgap2 > 0 & nrow(merged) != 1){
      tbmerged <- merged
      while (T){
        merged <- data.frame()
        i <- 1
        while (i<(nrow(tbmerged))){ # for loop which we can quickly skip through
          next_one <- which(tbmerged[(i+1):nrow(tbmerged),3]<=(tbmerged[i,4]+(tbmerged[i,4]-tbmerged[i,3])*maxgap2))
          if(length(next_one)>0){ # if true then that means we have overlap
            tmp <- which(tbmerged[,4]==max(c(tbmerged[next_one+i,4],tbmerged[i,4])))[1]
            merged <- rbind(merged,cbind(tbmerged[i,1],tbmerged[i,2],tbmerged[i,3],tbmerged[tmp,4]), stringsAsFactors = FALSE)
            i <- i+max(next_one)+1 # skip to the next not-overlapping sequence
            if (i==(nrow(tbmerged))){
              lastmerged <- FALSE
            } else {
              lastmerged <- TRUE
            }
          } else { # this sequence is not overlapping with others
            merged <- rbind(merged,cbind(tbmerged[i,1],tbmerged[i,2],tbmerged[i,3],tbmerged[i,4]), stringsAsFactors = FALSE)
            i <- i+1
            lastmerged <- FALSE
          }
        }
        if (!lastmerged){
          merged <- rbind(merged,cbind(tbmerged[nrow(tbmerged),1],tbmerged[nrow(tbmerged),2],tbmerged[nrow(tbmerged),3],tbmerged[nrow(tbmerged),4]), stringsAsFactors = FALSE)
        }
        merged <- as.data.frame(merged, stringsAsFactors = F)
        merged[,3] <- as.numeric(merged[,3])
        merged[,4] <- as.numeric(merged[,4])
        if(nrow(tbmerged)==nrow(merged) | nrow(merged) == 1){break}
        tbmerged <- merged
      }
    }
    
    #print("merge step 3")
    if (datatype == 1){
      colnames(merged) <- c("ref_name","query_name","query_start","query_end")
    } else if (datatype == 2){
      colnames(merged) <- c("ref_name","query_name","ref_start","ref_end")
    } else if (datatype == 3){
      colnames(merged) <- c("ref_name","query_name","query_start","query_end","ref_start","ref_end",colnames(data)[c(6,5,7)]) # get everything but still use query formatting
    } else if (datatype == 4){
      colnames(merged) <- c("ref_name","query_name","ref_start","ref_end","query_start","query_end",colnames(data)[c(5,6,7)]) # get everything but still use ref formatting
    } else if (datatype == 5){
      colnames(merged) <- colnames(data[,-5])
    } else {
      stop("datatype needs to be either 1,2,3,4, or 5")
    }
    
  } else {stop("Data must be either NA or a numeric value >= 0")}
  if (mode == 0){
    write.table(merged, gsub(".tsv",".NR.tsv",filename), row.names = F, sep = "\t")
  } else {
    return(merged)
  }
}

# calculate the percentage of aligned bases
# if merge gets replaced by bed tools then this would best be performed in the same script as merge, which would best be made in python
get_aligned_perc <- function(filename = "Temp/filteredSuper-Scaffold_1000002-H12.NR.tsv", fastafile = "Super-Scaffold_1000001.fasta", mode = 1, data = merged, out = opt$out){
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
    SS_seq <- readDNAStringSet(fastafile)
    
    #print(data)
    aligned_perc <- cbind(ref_name = data$ref_name[1],
                          query_name = data$query_name[1],
                          aligned_perc = sum(data$length)/(nchar(SS_seq))*100) # length divided by total SS length *100 is % aligned bases
    #print(aligned_perc)
  } else {stop("Data must be either NA or a numeric value >= 0")}
  write.table(aligned_perc, out, sep = "\t", row.names = FALSE)
}
