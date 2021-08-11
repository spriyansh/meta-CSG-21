# Loading the library
library(phyloseq)
library(tidyverse)
library(ape)

######################
### Making Object ####
######################

# Reading taxonomic annotations
otu_tax <- read.csv("asv_silva_tax_raw.tsv", sep = "\t")

# Column to rownames
otu_tax <- otu_tax %>% tibble::column_to_rownames("ASV")

# Conversion to matrix
otu_tax <-  as.matrix(otu_tax)

# Reading Sequence Variants
otu_abundance <- read.csv("ASV_abundance_raw.tsv", sep = "\t")

# Column to rownames
otu_abundance <- otu_abundance %>% tibble::column_to_rownames("ASV")

# Conversion to matrix
otu_abundance <-  as.matrix(otu_abundance)

# Reading Metadata
metadata <- read.csv("Meta_data.tsv", sep = "\t")

# Adding X to s_number
metadata$s_number <- sub("^", "X", metadata$s_number)

# Column to rownames
metadata  <- metadata %>%  tibble::column_to_rownames("s_number") 

# Reading the tree
tree <- read.tree("asv_tree.tree")

# Individual Objects
otu <- otu_table(otu_abundance, taxa_are_rows = T)
Taxa <- tax_table(otu_tax)
Samples <- sample_data(metadata)

# Creating Object
otu_PS <- phyloseq(otu, Taxa, Samples,phy_tree(tree))

# Summary
otu_PS

# Saving to RDS
saveRDS(otu_PS, "ASV_PS_Raw.RDS")

###########################
#### PhyloSeq Analysis ####
###########################

# Viewing the Available ranks
rank_names(otu_PS) 

# Creating a table for number of features for each phyla
table(tax_table(otu_PS)[, "Phylum"], exclude = NULL)

# Removing the NAs and uncharacterized Phylums
temp_df <- as.data.frame(tax_table(otu_PS))
unclassified <- subset(temp_df, grepl("unclassified", temp_df$Phylum))$Phylum

ps0 <- subset_taxa(otu_PS, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Summary Before and After
otu_PS
ps0

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

# Viewing
head(prevdf)

# Computing average and total prevalence
temp_df <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
temp_df

# Vector to be removed
filterPhyla <- temp_df[temp_df$`2` < 100,]$Phylum

# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)

# Summary
otu_PS
ps0
ps1

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

## [1] 3.3 = 3

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)

# Summary
otu_PS
ps0
ps1
ps2

saveRDS(ps2, "ASV_PS2.RDS")

# How many genera would be present after filtering?
length(get_taxa_unique(ps1, taxonomic.rank = "Genus"))

ps3 = tax_glom(ps2, "Genus", NArm = TRUE)

# Summary
otu_PS
ps0
ps1
ps2
ps3

saveRDS(ps3, "ASV_PS3.RDS")
