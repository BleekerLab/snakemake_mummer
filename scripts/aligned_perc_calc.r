# Wrapper for the aligned_perc_calc_functions.r script

# Load packages
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser(description='Calculate percentage of aligned bases')
parser$add_argument("-f", "--filename", metavar = "character", type = "character", nargs="+", help = "list of input files")
parser$add_argument("--fasta", metavar = "character", type = "character", nargs="+", help = "list of query fasta files")
parser$add_argument("--out", metavar = "character", type = "character", nargs="+", help = "list of output files")
parser$add_argument("-l", "--length_threshold", metavar = "integer", type = "integer", default =  1000, help = "threshold value for length filter [default]")
parser$add_argument("-i", "--identity_threshold", metavar = "double", type = "double", default = 90.0, help = "threshold value for identity filter [default]")

opt <- parser$parse_args()

# Load functions
source("scripts/aligned_perc_calc_functions.r")

mode <- 1 # mode 1 is wrapper mode, 0 is standalone mode
# remake of aligned_perc_calc.py
oldtime <- Sys.time() # test speed of script
coords_files <- opt$filename
for (i in coords_files){
  filtered_data <- get_filtered_data(filename = i,
                                     length_threshold = opt$length_threshold,
                                     identity_threshold = opt$identity_threshold,
                                     mode = mode) # filter the data

  merged <- merge_data(filename = gsub(".coords",".tsv",gsub("sorted","Temp/filtered",i)),mode = mode, data = filtered_data) # remove redundant sequences and merge overlapping sequences
  get_aligned_perc(filename = gsub(".coords",".NR.tsv",gsub("sorted","Temp/filtered",i)),mode = mode, data = merged, fastafile = opt$fasta, out = opt$out) # calculate the percentage of aligned bases
}
#print(Sys.time()-oldtime)
