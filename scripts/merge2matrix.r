# Final step in getting the matrix of aligned bases

#print("does it even start?")
# Load packages
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
  tmp_colname <- aligned_perc_matrix[1,1]

  #print(aligned_perc_matrix)
  aligned_perc_matrix <- spread(as.data.frame(aligned_perc_matrix),"ref_name","aligned_perc") # spread out the data based on the chr number
  tmp_names <- aligned_perc_matrix[,1] # set column containing names as rownames
  aligned_perc_matrix <- aligned_perc_matrix[,-1] # remove column containing names
  if(length(data)>1){
    aligned_perc_matrix <- aligned_perc_matrix[,order(colnames(aligned_perc_matrix))] # order the columns
  } else {
    aligned_perc_matrix <- as.data.frame(aligned_perc_matrix, row.names = tmp_names, stringsAsFactors = False)
    colnames(aligned_perc_matrix) <- tmp_colname
  }
  write.table(aligned_perc_matrix,outfile,sep = "\t",quote=F) # write as .tsv for use in future scripts
  #return(aligned_perc_matrix)
}

merge2matrix(data = opt$filename, outfile = opt$out)
