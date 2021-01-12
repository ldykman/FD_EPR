### FUNCTIONAL DISSIMILARITY AND DENDROGRAM ###
# Analysis for "Functional Traits Provide New Insight into Recovery and Succession at Deep-sea Hydrothermal Vents"
# Lauren Dykman
# August, 2020

# Creates a dendrogram for using a csv data table of trait modalities, with rows being species and columns being traits.
# Calculates Gower Distance.
# Can handle numerical, ordinal, or categorical data.
# Changes data type of four columns to ordered categories.
# Can handle missing values.

rm(list = ls())

### INSTALLING PACKAGES ###

#install.packages("FD")
#install.packages("cluster")
#install.packages("dplyr")

library("FD")
library("cluster")
library("dplyr")

### SETTING WORKING DIRECTORY ###

path <- "/Users/laurendykman/Desktop/github/FD_EPR"
path_out <- "/Users/laurendykman/Desktop/github/FD_EPR/output"

getwd()
setwd(path)
getwd()

### IMPORTING TRAIT TABLE ###

traits.raw <- read.csv("DataS2_BCO-DMO_Dykman.csv", row.names = 1, check.names = FALSE) # Importing trait data
traits.raw <- within(traits.raw, rm("FUNCTIONAL_GUILD")) # Removing functional group column for the dissimilarity calculation
traits <- dplyr::select(traits.raw, !dplyr::contains("CITATION")) # Removing citation columns
traits <- dplyr::select(traits, !dplyr::contains("aphiaID")) # Removing TAXON APHIA ID column
traits <- dplyr::select(traits, !dplyr::contains("scientificNameID")) # Removing TAXON APHIA ID column
str(traits) # Checking the orignial column type

### ADJUSTING COLUMNS TYPE TO NUMERICAL, CATEGORICAL, OR ORDERED ###

traits$MAXIMUM_ADULT_BODY_SIZE <- ordered(traits$MAXIMUM_ADULT_BODY_SIZE, levels=c('Small (~1mm)','Medium (~10mm)','Large (~100mm)','Very large (~1000mm)'))
traits$TROPHIC_MODE <- ordered(traits$TROPHIC_MODE, levels=c('Symbiont', 'Bacterivore', 'Detritivore', 'Carnivore S', 'Carnivore O'))
traits$RELATIVE_ADULT_MOBILITY <- ordered(traits$RELATIVE_ADULT_MOBILITY, levels=c('Sessile', 'Movement restricted', 'Crawler', 'Freely mobile'))
traits$EXTERNAL_PROTECTION <- ordered(traits$EXTERNAL_PROTECTION, levels=c('Soft bodied', 'Moderately protected', 'Well protected'))
traits[traits=="Unknown"] <- NA # Replace Unknown with NA

# IMPORTANT! Make sure all trait categories have changed as desired, and all cells still contain values.

str(traits)

# MAXIMUM ADULT BODY SIZE should be type ordered
# TROPHIC MODE should be type ordered
# RELATIVE ADULT MOBILITY should be type ordered
# EXTERNAL PROTECTION should be type ordered

### GENERATING GOWER DISSIMILARITY FROM THE TRAIT TABLE ###
# Gower dissimilarity is used because it can take numerical, categorical, or ordered data types. It can also accept missing data (Legendere and Legendere, 1998)

gow <- as.dist(gowdis(traits, ord = c("podani")))
write.csv(as.matrix(gow), file=paste(path_out, 'FD_EPR_GowerDis.csv', sep='/')) # Un-comment to save gower dissimilarity matrix as a .csv file

# CREATING A DENDROGRAM

cluster_dist <- hclust(gow, method="ward.D2") # Generating clusters using method Ward.D2 which implements Ward's (1963) clustering criterion (Murtagh and Legendre 2014)
dendro <- as.dendrogram(cluster_dist)
pdf(paste(path_out, paste("FD_EPR_DENDROGRAM_",Sys.Date(),".pdf"), sep='/'), height=7, width=11) # Saving the dendrogram as a figure
par(mar=c(mar=c(20,5,4,1)))
par(cex=0.9)
plot(dendro, ylab="Dissimilarity")
rect.hclust(cluster_dist, h=0.22, border='blue')
par(cex=1)
dev.off()

# OPTIMIZING CLUSTERS

clust.means <- c()
for (n in 2:40) {
  clust <-  cutree(cluster_dist, k=n)
  si <- silhouette(clust, gow)
  clust.means <- c(clust.means, mean(si[,3]))
}
plot(2:40, clust.means)
greatest <- which.max(clust.means) + 1
optimum <- silhouette(cutree(cluster_dist, k=greatest), gow)

# SAVING FUNCTIONAL GUILDS AS NEW COLUMN IN TRAIT TABLE

guilds <- optimum[,1]
letter.options <- LETTERS[1:length(guilds)]
letter.guilds <- c()
for (i in guilds) {
  letter.guilds <- c(letter.guilds, letter.options[i])
}
traits.new <- cbind(traits.raw, letter.guilds)
names(traits.new)[as.numeric(length(names(traits.new)))] <- "FUNCTIONAL_GUILD"
write.csv(as.matrix(traits.new), file=paste("DataS2_BCO-DMO_Dykman.csv"))

