# Get the locations of N nucleotide sequences in a DNA sequence

# Load packages
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description='Get the locations of N nucleotide sequences in a DNA sequence')
parser$add_argument("--fasta", metavar = "character", type = "character", help = "fasta file")
parser$add_argument("--out", metavar = "character", type = "character", nargs="+", help = "list of output files")
parser$add_argument("--n_threshold", metavar = "integer", type = "integer", default = 0, help = "threshold value for N filter [default]")

opt <- parser$parse_args()

# Get the locations on the DNA sequence that contain Ns
get_N_locations <- function(fastafile = opt$fasta, N_threshold = opt$n_threshold){
  my_seq <- readDNAStringSet(fastafile) 
  N_locations <- unlist(gregexpr("N",my_seq))
  if (N_locations[1] < 0){
    list_of_n_sequences <- cbind(start = NA, end = NA, length = NA)
  } else {
    sequence_mode <- F # is true while we're looking for a sequence of Ns instead of a single N
    list_of_n_sequences <- c()
    for (i in 1:(length(N_locations)-1)){
      if(N_locations[i]+1==N_locations[i+1]){
        if(!sequence_mode){
          mystart <- N_locations[i]
          sequence_mode <- T
        }
      } else {
        if(sequence_mode){
          myend <- N_locations[i]
          list_of_n_sequences <- rbind(list_of_n_sequences,c(start = mystart, end = myend)) # found a whole sequence of Ns
          sequence_mode <- F
        } else {
          list_of_n_sequences <- rbind(list_of_n_sequences,c(start = N_locations[i], end = N_locations[i])) # found a single N
          sequence_mode <- F
        }
      } 
    }
    if(sequence_mode){
      list_of_n_sequences <- rbind(list_of_n_sequences,c(start = mystart, end = N_locations[i+1]))
    }
    list_of_n_sequences <- cbind(list_of_n_sequences, length = list_of_n_sequences[,2]-list_of_n_sequences[,1]+1)
    if (N_threshold>0){
      list_of_n_sequences <- list_of_n_sequences[-which(list_of_n_sequences[,3]<N_threshold),]
    }
    if (length(dim(list_of_n_sequences))!=2){
      list_of_n_sequences <- cbind(start = list_of_n_sequences[1], end = list_of_n_sequences[2], length = list_of_n_sequences[3])
    }
  }
  return(list_of_n_sequences)
}

tmp <- get_N_locations(fastafile = opt$fasta)
write.table(tmp,file = opt$out, sep = "\t")

