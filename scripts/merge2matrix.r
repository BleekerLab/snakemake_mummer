# Final step in getting the matrix of aligned bases

# Load packages
#.libPaths("/home/wotten/R/x86_64-pc-linux-gnu-library/3.4")
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(argparse))

# get command line input
parser <- ArgumentParser(description='create final results.tsv')
parser$add_argument("-f", "--filename", metavar = "character", type = "character", nargs="+", help = "list of input files")
parser$add_argument("-o", "--out", metavar = "character", type = "character", default =  "results.tsv", help = "filename of output [default]")
opt <- parser$parse_args()


# Merge the aligned_percentages of all alignments into 1 matrix and writes it as .tsv
merge2matrix <- function(data = opt$filename, outfile = opt$out){
  aligned_perc_matrix <- c()
  for (i in data){
    aligned_perc_matrix <- rbind(aligned_perc_matrix,read.table(i, header = T, sep = "\t", stringsAsFactors = F))
  }
#  print(aligned_perc_matrix)
  aligned_perc_matrix <- spread(as.data.frame(aligned_perc_matrix),"ref_name","aligned_perc") # spread out the data based on the chr number
  rownames(aligned_perc_matrix) <- aligned_perc_matrix[,1] # set column containing names as rownames
  aligned_perc_matrix <- aligned_perc_matrix[,-1] # remove column containing names
  aligned_perc_matrix <- aligned_perc_matrix[,order(colnames(aligned_perc_matrix))] # order the columns
  write.table(aligned_perc_matrix,outfile,sep = "\t") # write as .tsv for use in future scripts
  #return(aligned_perc_matrix)
}

merge2matrix(data = opt$filename, outfile = opt$out)
