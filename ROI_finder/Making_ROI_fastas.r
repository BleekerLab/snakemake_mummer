# Make sure you have the correct file system setup.
# Most files belong in the input folder


# Source this file to load required packages and custom functions
source("ROI_fasta_functions.r") # you can generally ignore the warnings

# set which chromosome the ROI is from, where its first and last basepair are and the desired name..
Chromosome <- 1
LowerBound <- 2588520
UpperBound <- 2738521
ROI_name <- "P450"


# get alignments of all super-scaffolds to the chromosome that contains your region of interest.
tmp_filelist <- list.files("Input/Filtered/", paste0("chromosome.",Chromosome,".filtered"))

# for best results, take only the super-scaffolds that were assigned to that chromosome.
tmp_table <- assign_SStoChr()
matched_scaffolds <- rownames(tmp_table)[which(tmp_table$Outliers == Chromosome)]
matched_scaffolds <- substring(matched_scaffolds,nchar("Super-Scaffold_")+1)
tmp_filelist <- tmp_filelist[grep(paste(matched_scaffolds,collapse = "|"),tmp_filelist)]

# find out which super-scaffold best matches the region.
total_matching_length <- c()
for (file in tmp_filelist){
  tmp_filtered <- read.table(paste0("Input/Filtered/",file), header = T)
  tmp_filtered <- tmp_filtered[order(tmp_filtered[,1]),]
  
  # only pick the alignments which overlap with the region of interest
  tmp_filtered <- tmp_filtered[which((tmp_filtered[,1] > LowerBound & tmp_filtered[,1] < UpperBound) | (tmp_filtered[,2] > LowerBound & tmp_filtered[,2] < UpperBound)),]
  
  # make a metric of how well each super-scaffold matches, in this case the sum of matchin alignment lengths
  total_matching_length <- c(total_matching_length, sum(tmp_filtered$length_1))
  
}
tmp_SS <- tmp_filelist[which(total_matching_length==max(total_matching_length))]
best_matching_SS <- substr(tmp_SS, 0, regexpr("_vs_",tmp_SS)-1)




# read the fasta of your reference
Heinz_ROI <- readDNAStringSet(paste0("Input/Ref_fastas/Solanum_lycopersicum.SL3.0.dna.chromosome.",Chromosome,".fasta"))
# select the part that you know the ROI to be in
Heinz_ROI <- substr(Heinz_ROI,LowerBound,UpperBound)

# write the ROI of the reference
write.fasta(Heinz_ROI,paste0("Heinz_chr",Chromosome,"_",LowerBound,"-",UpperBound,"_",ROI_name),paste0("Output/Heinz_",ROI_name,".fasta"))


# read the filtered file of the best fitting alignment to that region
# I just know that this is the right super scaffold for this ROI but for others you need to perform an additional check not included here
filtered_alignment <- read.table(paste0("Input/filtered/",best_matching_SS,"_vs_Solanum_lycopersicum.SL3.0.dna.chromosome.",Chromosome,".filtered"), sep = "\t", header = T)
tmp_table <- filtered_alignment

# Need to check if it is sorted or not
for (i in 1:nrow(tmp_table)){
  if (tmp_table[i,1]>tmp_table[i,2]){
    print("it is not sorted")
    break()
  }
}

# only get hits to our ROI
# the start of the hit OR the end of the hit needs to be within our ROI
filtered_alignment <- tmp_table[which((tmp_table[,1] > LowerBound & tmp_table[,1] < UpperBound) | (tmp_table[,2] > LowerBound & tmp_table[,2] < UpperBound)),]

# manually check which entries are part of the big cluster of hits, remove the rest.
# in this case you can make a histogram to clearly see 2 outliers
hist(filtered_alignment[,3],xlab = "alignment_location", main = "distribution of alignment locations")
View(filtered_alignment) # manually look up which hits are outside of the area with the highest density of hits

# If you have found a manual cut off, insert those here, if you think you don't need any, set both to 0..
cutoff_upper <- 0
cutoff_lower <- 0
if (cutoff_upper==0 & cutoff_lower==0){
  ROI_hits <- filtered_alignment
} else {
  ROI_hits <- filtered_alignment[which(filtered_alignment[,3]< cutoff_upper & filtered_alignment[,4]>cutoff_lower),] # or use a cutoff value based on the histogram
}

# extract the range where you expect to find the ROI
ROI_bounds <- c(min(ROI_hits[,3:4]),max(ROI_hits[,3:4]))

# create fasta
SS_ROI <- readDNAStringSet(paste0("Input/Query_fastas/",best_matching_SS,".fasta"))
SS_ROI <- substr(SS_ROI,ROI_bounds[1],ROI_bounds[2])
write.fasta(SS_ROI,paste0(best_matching_SS,"_", ROI_bounds[1],"-", ROI_bounds[2],"_",ROI_name),paste0("Output/",best_matching_SS,"_LTP.fasta"))

# OPTIONAL: check if the range is as expected
(UpperBound-LowerBound)/(ROI_bounds[2]-ROI_bounds[1])
# If this value is a lot lower than 1 then there might be a large insert of N in the query.
# In that case extention corrections won't work.

# If this value is greater than 1 then it'd be useful to extend the ROI_bounds for the super-scaffold.
# Base the extention on the previously calculated value.
if ((UpperBound-LowerBound)/(ROI_bounds[2]-ROI_bounds[1])>1){
  extention <- (UpperBound-LowerBound)/(ROI_bounds[2]-ROI_bounds[1])-1
  ROI_bounds <- c(ROI_bounds[1]-(ROI_bounds[2]-ROI_bounds[1])*extention, ROI_bounds[2]+(ROI_bounds[2]-ROI_bounds[1])*extention)
}

# if we check the range again we can see that it is much closer to expected
(UpperBound-LowerBound)/(ROI_bounds[2]-ROI_bounds[1])

# and correct the earlier created fasta
SS_ROI <- readDNAStringSet(paste0("Input/Query_fastas/",best_matching_SS,".fasta"))
SS_ROI <- substr(SS_ROI,ROI_bounds[1],ROI_bounds[2])
write.fasta(SS_ROI,paste0(best_matching_SS,"_", ROI_bounds[1],"-", ROI_bounds[2],"_LTP_cluster"),paste0("Output/",best_matching_SS,"_LTP.fasta"))


