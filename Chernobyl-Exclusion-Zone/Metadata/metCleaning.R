##################################################################
###### Cleaning metadata, renaming and Zipping fastq files #######
##################################################################

##################################################################
# Author Information #############################################
# Name # Priyansh Srivastava #####################################
# Contact # spriyansh29@gmail.com ################################
##################################################################

# Imports
library(tidyverse)

# Loading the metadata
metadataRaw <- read.csv("md_for_phyloseq_Priyansh.csv", header = F)

# Viewing the metadata
#View(metadataRaw)

# Saving subset
metadataSubset <- metadataRaw[c(8:74), ]

# Renaming the columns
colnames(metadataSubset) <- metadataSubset[1,]
metadataSubset <- metadataSubset[-1, ] 

# Reordering by sample numbers
metadataSubset <- metadataSubset[order(as.numeric(metadataSubset$s_number)),]

# Adding custom names in column: neo_names
metadataSubset$neo_names <- paste(metadataSubset$s_number, metadataSubset$s_name, 
                                  sep = "_g")
metadataSubset$neo_names <- gsub("(\\_.*?)", "\\_r", metadataSubset$neo_names)
metadataSubset$neo_names <- str_replace(metadataSubset$neo_names, "rg", "g")

# Reset index
row.names(metadataSubset) <- NULL

# Viewing the subset
View(metadataSubset)

###############################################
#### Designing New Dataframe for merging ######
###############################################

# Reading all the filenames
allFileNamesOld <- list.files(pattern = "*001.fastq*")

# Removing patterns "L001_" and 
allFileNamesTemp <- str_remove_all(allFileNamesOld, "L001_")
allFileNames <- str_remove_all(allFileNamesTemp, "_001")

# Replacing patterns "-" with "g_" only first occurrence
allFileNamesTemp2 <- str_replace(allFileNames, "-", "_g")

# Replacing the "." to run
allFileNames <- paste0(sub('\\.', '_r', allFileNamesTemp2))

# Correcting for blank and Mock and .fastq
allFileNames <- str_replace(allFileNames, "R1_r", "R1.") # for R1
allFileNames <- str_replace(allFileNames, "R2_r", "R2.") # for R2
allFileNames <- str_replace(allFileNames, "gB", "gBlank_r") # For B
allFileNames <- str_replace(allFileNames, "Mock_Community", "Mock") # for Mock

# Regular Expressions for subsetting adaptors
regObjectAdaptors <- regexpr("[A|G|T|C](.*)[A|G|T|C]", allFileNames)
adaptors <- regmatches(allFileNames, regObjectAdaptors)

# Regular Expressions for pairs
regObjectPairs <- regexpr("R(.)", allFileNames)
pairs <- regmatches(allFileNames, regObjectPairs)

# New File names
names <- str_replace(allFileNames, "_[A|G|T|C](.*)[A|G|T|C]_", "_")

# Merge Column
regObjectNeoMerge <- regexpr("(*.*)_[A|G|T|C]", allFileNames)
NeoMerge <- regmatches(allFileNames, regObjectNeoMerge)
NeoMerge <- str_replace(NeoMerge, "_[A|G|T|C]", "")

# Making the dataframe for old and new filenames 
neoMetadata <- data.frame(NeoMerge, allFileNamesOld, names, adaptors, pairs)

# Renaming the column
colnames(neoMetadata) <- c("neo_names","Old_File_Names", "New_file_Names",
                           "Adaptors", "Pair")

# Viewing the new metadata
#View(neoMetadata)

# Merging the old and new metadata
metadata <- merge(metadataSubset, neoMetadata, by="neo_names")

# Viewing
View(metadata)

# Writing
write.csv(metadata, "finalMetadata.csv")

###############################################
################ Missing data #################
# Missing for sample number 5, 31 and 57 ######
# Missing metatdata for group 2, 11 and 24 ####
# Missing for Mock and black files ############
# Solution Mannually appending ################
# Out of 146 files 14 files have no metadata ##
# Proceeding with 132 files ################## 
###############################################

