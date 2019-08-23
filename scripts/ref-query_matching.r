# merging and ref-query matching script for the smakemake mummer pipeline (by Bleeker group)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser(description='Calculate percentage of aligned bases')
parser$add_argument("-f", "--filename", metavar = "character", type = "character", help = "list of input files")
parser$add_argument("--fasta", metavar = "character", type = "character", nargs="+", help = "list of query fasta files")
parser$add_argument("--out", metavar = "character", type = "character", nargs="+", help = "list of output files")

opt <- parser$parse_args()

# Load functions
source("scripts/aligned_perc_calc_functions.r")

N_filtered_data = read.table(opt$filename, sep = "\t", header = T, stringsAsFactors = FALSE)
merged <- merge_data(data = N_filtered_data) # remove redundant sequences and merge overlapping sequences
get_aligned_perc(data = merged, fastafile = opt$fasta, out = opt$out) # calculate the percentage of aligned bases
