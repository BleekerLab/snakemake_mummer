# First step in getting the matrix of aligned bases. Merging inputs from the queries

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(argparse))

# get command line input
parser <- ArgumentParser(description='create final results.tsv')
parser$add_argument("-f", "--filename", metavar = "character", type = "character", nargs="+", help = "list of input files")
parser$add_argument("-o", "--out", metavar = "character", type = "character", default =  "results.tsv", help = "filename of output [default]")
opt <- parser$parse_args()


aligned_perc_matrix <- c()
for (i in opt$filename){
  aligned_perc_matrix <- rbind(aligned_perc_matrix,read.table(i, header = T, sep = "\t", stringsAsFactors = F))
}
write.table(aligned_perc_matrix,opt$out,sep = "\t",quote=F)