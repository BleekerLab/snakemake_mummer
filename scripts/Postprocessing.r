# Load packages
library(Biostrings)
library(ggplot2)
library(seqinr)

# Load functions
source("scripts/postprocessing_functions.r")

# Loading the file
aligned_perc_matrix <- get_aligned_perc_matrix(filename = "fulldataresultsV0.1.1a.tsv")
mylist <- get_chr_identifiers(aligned_perc_matrix = aligned_perc_matrix)
chrID <- mylist$chrID
commonname <- mylist$commonname
colnames(aligned_perc_matrix) <- chrID
aligned_perc_matrix <- aligned_perc_matrix[,order(as.numeric(colnames(aligned_perc_matrix)))]

aligned_perc_matrix <- cbind(aligned_perc_matrix, best = apply(aligned_perc_matrix,1,max))
for (i in 1:nrow(aligned_perc_matrix)){
  aligned_perc_matrix[i,13] <- which(aligned_perc_matrix[i,-13] == aligned_perc_matrix[i,13])
}
coords_files <- paste0("scratch/coords/", rownames(aligned_perc_matrix[which(aligned_perc_matrix[,13]==1),]), "_vs_", commonname, 1,".sorted.coords")


coords_files <- list.files("scratch/coords/")
filename <- paste0("scratch/coords/",coords_files[1])


source("scripts/aligned_perc_calc_functions.r")


# filter and merg all sorted coords entries for SS that best matched a certain chromosome
opt = list(length_threshold = 1000, identity_threshold = 85.0, fasta = "SSvsHeinz/queries/Super-Scaffold_1000033.fasta", n_threshold = 10, out = "testfile.tsv")
chr_1_alignments <- c()
for (filename in coords_files){
  filtered_data <- get_filtered_data(filename = filename,
                                     length_threshold = opt$length_threshold,
                                     identity_threshold = opt$identity_threshold,
                                     mode = 1) # filter the data
  N_file <- paste0("scratch/SSvsHeinz/query_N/", substr(filename,gregexpr("coords/",filename)[[1]][1]+nchar("coords/"),gregexpr("_vs_",filename)[[1]][1]-1), ".txt")
  N_filtered_data <- filter_N_entries(data = filtered_data, N_file = N_file)
  chr_1_alignments <- rbind(chr_1_alignments, N_filtered_data)
}
chr_1_alignments <- chr_1_alignments[order(chr_1_alignments[,1]),]
chr_1_merged <- merge_data(filename = gsub(".coords",".tsv",gsub("sorted","Temp/filtered",i)),mode = 1, data = chr_1_alignments, datatype = 2)
chr_1_merged <- cbind(chr_1_merged,chr_1_merged[,4]-chr_1_merged[,3])
chr_seq <- readDNAStringSet("SSvsHeinz/refs/Solanum_lycopersicum.SL3.0.dna.chromosome.1.fasta")
sum(chr_1_merged[,5])/(nchar(chr_seq))


filter_N_entries <- function(data = filtered_data, N_file = opt$nfile){
  tmp <- read.table(N_file, sep="\t")
  N_containing_rows <- c()
  N_Ndata <- c()
  for (i in 1:nrow(tmp)){
    if (any(data[,3]<tmp[i,1] & data[,4]>tmp[i,1])){
      N_containing_rows <- c(N_containing_rows,which(data[,3]<tmp[i,1] & data[,4]>tmp[i,1]))
      N_Ndata <- c(N_Ndata, i)
    }
  }
  N_filtered_data <- data[-N_containing_rows,]
  return(N_filtered_data)
}

# Making plot with the amount of SS aligned per chr
# best_matches <- c()
# for (i in 1:nrow(aligned_perc_matrix)){
#  best_matches <- c(best_matches,colnames(aligned_perc_matrix)[which(aligned_perc_matrix[i,]==max(aligned_perc_matrix[i,]))])
# }
# names(best_matches) <- rownames(aligned_perc_matrix)
# for (i in 1:length(unique(best_matches))){
#   new_name <- substring(unique(best_matches)[i],gregexpr("[A-Z]",unique(best_matches)[i])[[1]]+1)
#   if (nchar(new_name)==1){
#     new_name <- paste0(0,new_name)
#   }
#   best_matches[which(best_matches==unique(best_matches)[i])] <- new_name
# }
# best_matches <- as.numeric(best_matches[order(best_matches)])
freq_best_matches <- c()
for (i in 1:max(aligned_perc_matrix[,13])){
  freq_best_matches <- c(freq_best_matches, sum(aligned_perc_matrix[,13]==i))
}
names(freq_best_matches) <- 1:12
barplot(freq_best_matches, ylab = "matches", xlab = "chromosome", main = "Super-Scaffolds aligned per chromsome")
for (i in 1:12){
  text(((i)*1.2)-0.5,2.5, labels = freq_best_matches[i])
}

# How much of each chromosome is actually aligned?
coords_files <- paste0("scratch/coords/",list.files("scratch/coords/"))
for (i in coords_files){
  opt = list(length_threshold = 1000, identity_threshold = 85.0, fasta = "SSvsHeinz/queries/Super-Scaffold_1000001.fasta", out = "testfile.tsv")
  filtered_data <- get_filtered_data(filename = i,
                                     length_threshold = opt$length_threshold,
                                     identity_threshold = opt$identity_threshold,
                                     mode = 1) # filter the data
  #N_filtered_data <- filter_N_entries(data = filtered_data, N_file = opt$out)
  
}




# Merge N but for actual alignments
coords <- read.table(filename, sep = "", skip = 5)
merged <- merge_data(mode = 1, data = filtered_data, maxgap = 2000)



# Load the chromosome
chr_seq <- readDNAStringSet(paste0("~/Wicher_comparitive_genomics/02SS_alignments_nucmer/Input/Heinz_chromosomes/Solanum_lycopersicum.SL3.0.dna.chromosome.",data$query_name[1],".fasta"))

hist(data$`%_identity`[which(data$length_1 >= 1000)])
hist(data$length_1[which(data$`%_identity` >= 99.0)])










# Whole fuckin lot of Ns
SS_fastafiles <- list.files("SSvsHeinz/queries/")
SS_N_locations <- get_N_locations(paste0("SSvsHeinz/queries/",SS_fastafiles))

chr_fastafiles <- list.files("~/Wicher_comparitive_genomics/02SS_alignments_nucmer/Input/Heinz_chromosomes/")
chr_N_locations <- get_N_locations(chr_fastafiles, datatype = 2)






# further analysis of N distribution on the DNA sequence
N_merged <- merge_Ns(data = SS_N_locations[[13]], maxgap = 5)
N_merged <- cbind(N_merged, length = N_merged[,2]-N_merged[,1]+1)

filename <- fastafiles[13]
my_seq <- readDNAStringSet(paste0("SSvsHeinz/queries/",filename))
filter_threshold <- 1
tmp <- N_merged[which(N_merged[,3]>=filter_threshold),]

myinterval <- round(nchar(my_seq)/10^(nchar(nchar(my_seq))-2))*10^(nchar(nchar(my_seq))-2)


plot.new()
axis(1, at = seq(0,1, 1/(myinterval/10^(nchar(nchar(my_seq))-1))), labels = round(seq(0,myinterval,10^(nchar(nchar(my_seq))-1))))
axis(2, at = seq(0,1, 0.2), labels = round(seq(0,max(tmp[,3]),max(tmp[,3])/5)))
xlocations <- tmp[,1:2]/myinterval
ylocations <- tmp[,3]/max(tmp[,3])
rect(xlocations[,1],0,xlocations[,2],ylocations)
mtext("position on the DNA sequence",1, line = 3)
mtext("length of N sequence",2, line = 3)
mtext(gsub(".fasta","", filename),3, line = 1, font = 2, cex = 1.25)

#plot(apply(tmp[,1:2],1,mean),tmp[,3], type = "h", xlab = "location on DNA sequence", ylab = "length of N sequence", main = gsub(".fasta","",filename))

#General data information
data_stats <- c()
queries <- list.files("SSvsHeinz/queries/")
for (i in queries){
  tmp_seq <- readDNAStringSet(paste0("SSvsHeinz/queries/", i))
  data_stats <- rbind(data_stats, cbind(DNA_seq = substr(i,1,regexpr(".fasta",i)[[1]][1]-1),
                                        nchar = nchar(tmp_seq),
                                        N = length(gregexpr("N",tmp_seq)[[1]])))
}
data_stats <- cbind(data_stats, N_percentage = round(100*as.numeric(data_stats[,3])/as.numeric(data_stats[,2]),2))




# # Actual postprocessing
# Out_align <- get_out_align() # get the names of the output of the aligned bases files
# Tab_align <- get_tab_align(filenames = Out_align)# Generate table with on row = chr and col = SS and values are % aligned bases
# 
# # Calculate the ratio between the highest and second highest aligned bases per SS
# tmp <- Tab_align
# for (i in 1:ncol(tmp)){
#   tmp[which(tmp[,i]==max(tmp[,i])),i] <- 0
# }
# best_2ndbest_ratios <- round((apply(Tab_align,2,max)-apply(tmp,2,max))/apply(Tab_align,2,max),2)
# 
# # Create histograms of max values and 2nd highest values
# par(mfrow=c(2,2))
# hist(apply(Tab_align,2,max),
#      xlim = c(0,100),
#      #ylim = c(0,20),
#      main = "best alignments")
# hist(c(tmp),
#      xlab = "2nd best",
#      xlim = c(0,100),
#      #ylim = c(0,20),
#      main = "non-best alignments")
# hist(apply(Tab_align,2,max)-apply(tmp,2,max),
#      xlab = "% alignment best - 2nd best",
#      xlim = c(0,100),
#      #ylim = c(0,20),
#      main = "best alignments - 2nd best")
# 
# # Create density plot
# best <- data.frame(value = apply(Tab_align,2,max))
# best2 <- data.frame(value = apply(tmp,2,max))
# best$kind <- "best"
# best2$kind <- "2nd best"
# tmptab <- rbind(best,best2)
# ggplot(tmptab, aes(value, fill = kind)) + geom_density(alpha = 0.2)
# 
# ggplot(tmptab, aes(value, fill = kind)) + 
#   geom_histogram(alpha = 0.5, aes(), position = 'identity')
# 
# 
