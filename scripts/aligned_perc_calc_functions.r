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
  
  # Read and fix
  data <- read.table(filename, sep = "", skip = 5)
  data <- data[,-which(data[1,]=="|")]
  colnames(data) <- tolower(c("Start_ref","End_ref","Start_query","End_query","Length_1","Length_2","%_Identity","Ref_name","Query_name"))
  
  # filter
  filtered_data <- data[-which(data$length_2 < length_threshold | data$`%_identity` < identity_threshold),]
  rownames(filtered_data) <- 1:nrow(filtered_data)
  if (mode == 0){
    write.table(filtered_data, paste0("Temp/filtered",filtered_data$ref_name[1],"-H",filtered_data$query_name[1],".tsv"), row.names = F, sep = "\t")
  } else if (mode == 1){
    return(filtered_data)
  }
  
}

# merge data
merge_data <- function(filename = "Temp/filteredSuper-Scaffold_1000002-H12.tsv",mode = 0, data = filtered_data){
  if (mode == 0){
    data <- read.table(filename, sep = "\t", header = T)
  }
  tbmerged <- data[,c(8,9,3,4)] # get just the relevant data
  tbmerged[,1] <- as.character(tbmerged[,1])
  
  while (T){
    merged <- c()
    i <- 1
    while (i<(nrow(tbmerged))){ # for loop which we can quickly skip through
      next_one <- which(tbmerged[(i+1):nrow(tbmerged),3]<=tbmerged[i,4])
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
  if (mode == 0){
    write.table(merged, gsub(".tsv",".NR.tsv",filename), row.names = F, sep = "\t")
  } else {
    return(merged)
  }
}

# calculate the percentage of aligned bases
get_aligned_perc <- function(filename = "Temp/filteredSuper-Scaffold_1000002-H12.NR.tsv", mode = 0, data = merged){
  if (mode == 0){
    data <- read.table(filename, sep = "\t", header = T, stringsAsFactors = F)
  }
  data <- cbind(data,length = (data[,4]-data[,3])) # calculate length
  SS_seq <- readDNAStringSet("test/queries/",paste0(data$query_name[1],".fasta")) 
  aligned_perc <- cbind(ref_name = data$ref_name[1],
                        query_name = data$query_name[1],
                        aligned_perc = sum(data$length)/nchar(SS_seq)*100) # length divided by total SS length *100 is % aligned bases
  write.table(aligned_perc, paste0(aligned_perc[,2],"_vs_",aligned_perc[,1],".txt"), sep = "\t", row.names = FALSE)
}
