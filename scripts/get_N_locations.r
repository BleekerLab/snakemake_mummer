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
  }
  return(longasslist)
}

tmp <- get_N_locations(fastafile = opt$fasta)
write.table(tmp,file = opt$out, sep = "\t")

