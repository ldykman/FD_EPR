### CALCULATING RAO'S QUADRATIC ENTROPY ###
# Lauren Dykman
# June 26, 2018

### INSTALLING PACKAGES ###

rm(list = ls())

#install.packages(c("FD", "dplyr", "vegan"))

library("FD")
library("dplyr")
library("vegan")

path <- "/Users/laurendykman/Desktop/github/FD_EPR"
path_out <- "/Users/laurendykman/Desktop/github/FD_EPR/output"

getwd()
setwd(path)
getwd()

### IMPORTING DATA ###

### IMPORTING TRAIT TABLE ###

traits.raw <- read.csv("DataS2_BCO-DMO_Dykman.csv", row.names = 1, check.names = FALSE) # Importing trait data
traits.raw <- within(traits.raw, rm("FUNCTIONAL_GUILD")) # Removing functional group column for the dissimilarity calculation
traits <- dplyr::select(traits.raw, !dplyr::contains("CITATION")) # Removing citation columns
traits <- dplyr::select(traits, !dplyr::contains("AphiaID")) # Removing TAXON APHIA ID column
traits <- dplyr::select(traits, !dplyr::contains("scientificNameID")) # Removing scientific name ID column
str(traits) # Checking the original column type

### ADJUSTING COLUMNS TYPE TO NUMERICAL, CATEGORICAL, OR ORDERED ###

traits$MAXIMUM_ADULT_BODY_SIZE <- ordered(traits$MAXIMUM_ADULT_BODY_SIZE, levels=c('Small (~1mm)','Medium (~10mm)','Large (~100mm)','Very large (~1000mm)'))
traits$HABITAT_COMPLEXITY <- as.factor(traits$HABITAT_COMPLEXITY)
traits$TROPHIC_MODE <- ordered(traits$TROPHIC_MODE, levels=c('Symbiont', 'Bacterivore', 'Detritivore', 'Carnivore S', 'Carnivore O'))
traits$FEEDING_METHOD <- as.factor(traits$FEEDING_METHOD)
traits$RELATIVE_ADULT_MOBILITY <- ordered(traits$RELATIVE_ADULT_MOBILITY, levels=c('Sessile', 'Movement restricted', 'Crawler', 'Freely mobile'))
traits$EXTERNAL_PROTECTION <- ordered(traits$EXTERNAL_PROTECTION, levels=c('Soft bodied', 'Moderately protected', 'Well protected'))
traits$LARVAL_DEVELOPMENT <- as.factor(traits$LARVAL_DEVELOPMENT)
traits$REPRODUCTIVE_TYPE <- as.factor(traits$REPRODUCTIVE_TYPE)
traits[traits=="Unknown"] <- NA # Replace Unknown with NA

# IMPORTANT! Make sure all trait categories have changed as desired, and all cells still contain values.

str(traits)

# MAXIMUM ADULT BODY SIZE should be type ordered
# TROPHIC MODE should be type ordered
# RELATIVE ADULT MOBILITY should be type ordered
# EXTERNAL PROTECTION should be type ordered
# IMPORTING ABUNDANCE AND SITE INFORMATION

### IMPORTING ABUNDANCE DATA AND GETTING ENVIRONMENTAL DATA ###

abund <- as.data.frame(read.csv("DataS1.csv", row.names = 1, stringsAsFactors=FALSE))

months.character <- as.character(as.factor(abund["Months Post-Eruption",]))
temperatures <- as.numeric(abund["Temperature",])
species <- rownames(abund)

abund <- abund[!rownames(abund) %in% c("Months Post-Eruption", "Temperature"), -1]
abund <- as.data.frame(lapply(abund, as.numeric))

rownames(abund) <- species[-c(1,2)]
colnames(abund) <- months.character[-1]

abund <- t(rowsum(t(abund), group = colnames(abund), na.rm = T)) # Sum columns from same time point
abund <- subset(abund, select=c(unique(months.character[-1]))) # Re-ordering columns
months.numeric <- c(-10, as.numeric(colnames(abund))[-1])
abund <- as.data.frame(t(abund))

#abund <- abund[, colnames(abund) != "X"]
#abund[abund == "Pre-Eruption"] <- as.numeric
#abund <- abund[!rownames(abund) %in% c("Temperature", "Months Post-Eruption"),]
#colnames(abund) <- as.character(months.character)
#abund <- abund[!rownames(abund) %in% c("Months Post-Eruption"),]
#colnames(abund)[colnames(abund) == "NA"] <- as.numeric(-10.0)

#abund <- sapply(abund, as.numeric)
#abund <- t(rowsum(t(abund), group = colnames(abund), na.rm = T)) # Sum columns from same time point


### SUBSETTING MATRICES ###
# If sum for any species is zero, write that species name in the vector called 'del' and those species names will be removed from analysis.

del <- c(colnames(abund)[colSums(abund)==0])
print(del)
traits <- traits[! row.names(traits) %in% del, drop = F, ]
abund <- abund[, ! colnames(abund) %in% del, drop = F]

### CHECKING DIMENSIONS ###
# The row total in traits should equal the column number in abund. In other words, all species names should be the same.

rownames(traits) == colnames(abund)
dim(traits)
dim(abund)

### FUNCTIONAL DIVERSITY ###

functdiv <- dbFD(traits, abund,
                ord = c("podani"), asym.bin = NULL,
                corr = c("none"),
                calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
                scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "kmeans",
                km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100,
                km.crit = c("calinski", "ssi"), calc.CWM = FALSE,
                calc.FDiv = TRUE, dist.bin = 2,
                print.pco = FALSE, messages = TRUE)

### FUNCTIONAL DIVERSITY OVER TIME ###

# Plotting RaoQ

par(mfrow=c(1,1))
RaoQ <- functdiv$RaoQ
plot(RaoQ~as.numeric(months.numeric), main = "RaoQ Over Time", xlab="Weeks After Eruption", ylab = "RaoQ")
linmod <- lm(RaoQ[2:length(RaoQ)]~as.numeric(months.numeric)[2:length(months.numeric)]+as.numeric(months.numeric**2)[2:length(months.numeric**2)])
fitx <- c(1:140)
fity <- linmod$coefficients[1] + linmod$coefficients[2]*fitx + linmod$coefficients[3]*fitx**2
lines(fitx, fity)
abline(v=0)
summary(linmod)
summary(linmod)$r.squared
f <- summary(linmod)$fstatistic
p <- pf(f[1],f[2],f[3],lower.tail=F)
attributes(p) <- NULL
p

write.csv(functdiv$RaoQ, paste(path_out, "FD_EPR_RaoQ.csv", sep="/"))

# Plotting Functional Richness

par(mfrow=c(1,1))
FRic <- functdiv$FRic
plot(FRic~as.numeric(months.numeric), main = "Functional Richness Over Time", xlab="Weeks After Eruption", ylab = "FRic")
linmod <- lm(FRic[2:length(FRic)]~as.numeric(months.numeric)[2:length(months.numeric)]+as.numeric(months.numeric**2)[2:length(months.numeric**2)])
fitx <- c(1:140)
fity <- linmod$coefficients[1] + linmod$coefficients[2]*fitx + linmod$coefficients[3]*fitx**2
lines(fitx, fity)
abline(v=0)
summary(linmod)
summary(linmod)$r.squared
f <- summary(linmod)$fstatistic
p <- pf(f[1],f[2],f[3],lower.tail=F)
attributes(p) <- NULL
p

write.csv(functdiv$FRic, paste(path_out, "FD_EPR_FRic.csv", sep="/"))

