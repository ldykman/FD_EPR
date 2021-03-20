### MULTINOMIAL LOGISTIC REGRESSION ###

# Analysis for "Functional Traits Provide New Insight into Recovery and Succession at Deep-sea Hydrothermal Vents"
# Lauren Dykman
# August, 2020

# Uses multinomial logistic regression in the package "nnet" to calculate coefficients and deviance for the relative changes in modalities within each trait as a function of time.
# This model is quadratic to allow for non-linear trends. P-values are calculated by 1,000 randomizations, which involves repeatedly randomizing the columns while keeping observation times fixed.
# Significance was assessed as the proportion of randomizations for which the randomized deviance was less than the observed deviance.

### LOADING PACKAGES AND SETTING WORKING DIRECTORY ###

# install.packages('foreign', 'nnet', 'lmtest', 'plyr', 'dplyr')

library('foreign')
library('nnet')
library('lmtest')
library('plyr')
library('dplyr')

rm(list = ls())

path <- "/Users/laurendykman/Desktop/github/FD_EPR"
path_out <- "/Users/laurendykman/Desktop/github/FD_EPR/output"

getwd()
setwd(path)
getwd()

### IMPORTING DATA ###

abund <- as.data.frame(read.csv("DataS1.csv", row.names = 1))
traits <- as.data.frame(read.csv("DataS2_BCO-DMO_Dykman.csv", row.names = 1, check.names = FALSE))
modality <- as.data.frame(read.csv("FD_EPR_Modality.csv"))

### MODIFYING DATA TABLES ###

abund <- abund[, abund[rownames(abund) == "Months Post-Eruption",] != "Pre-Eruption" & colnames(abund) != "X"] # Removing pre-eruption samples
traits <- dplyr::select(traits, !dplyr::contains("CITATION")) # Removing citation columns
traits <- dplyr::select(traits, !dplyr::contains("aphiaID")) # Removing TAXON APHIA ID column
traits <- dplyr::select(traits, !dplyr::contains("scientificNameID")) # Removing scientific name ID
str(traits) # Checking the original column type

### EXTRACT MONTHS DATA AND REMOVE ENVIRONMENTAL DATA FROM DATA FRAME ###

months.character <- as.character(abund["Months Post-Eruption",])
abund <- abund[!rownames(abund) %in% c("Temperature", "Months Post-Eruption"),]
colnames(abund) <- months.character
abund <- t(rowsum(t(abund), group = colnames(abund), na.rm = T)) # Sum columns from same time point
abund <- as.data.frame(subset(abund, select=c(unique(months.character)))) # Re-ordering columns
months.numeric <- as.numeric(colnames(abund))
months.character <- as.character(colnames(abund))

### PREPARE STATISTICS DATA FRAME ###

stats.baseline <- c("Small (~1mm)", "Does not add complexity", "Symbiont", "Non-feeding", "Sessile", "Soft bodied", "Lecithotrophic", "Gonochoristic", "A")   
stats.categories <- as.data.frame(stats.baseline, row.names = colnames(traits))
colnames(stats.categories) <- c("BASELINE")

### CHECK SPECIES MATCH ###

species.abund <- row.names(abund)
species.traits <- row.names(traits)
identical(species.abund, species.traits)

### CREATING DATA FRAME IN MULTINOM FORMAT ###

epr.data <- data.frame()
for(month in colnames(abund)){
  multinom.format <- as.data.frame(lapply(traits, rep, abund[,month]), check.names = FALSE)
  multinom.format <- cbind(rep(as.numeric(month), dim(multinom.format)[1]), multinom.format)
  epr.data <- rbind(epr.data, multinom.format)
}
colnames(epr.data)[1] <- "MONTHS"

# CREATING EXPLANATORY COLUMNS

epr.data$MONTHS.SQUARED = as.numeric(epr.data$MONTHS**2)
epr.data$MONTHS.FACTOR = as.factor(epr.data$MONTHS)

# LOOPING THROUGH TRAITS

traits.predicted <- data.frame()
all.deviance <- c()
all.deviance.factor <- c()
all.df <- c()
all.chi.squared <- c()
all.p.values <- c()

for (trait in colnames(traits)){
  
  # CALCULATING DEVIANCE
  
  epr.data$X <- droplevels(relevel(as.factor(epr.data[,trait]), ref = as.character(stats.categories[trait,1])))
  multinom.observed <- multinom(X~MONTHS+MONTHS.SQUARED, data = epr.data, maxit=800)
  multinom.null <- multinom(X~1, data = epr.data, maxit=800)
  multinom.factor <- multinom(X~MONTHS.FACTOR, data = epr.data, maxit=800)
  
  # TESTING SIGNIFICANCE
  
  linear.test <- lrtest(multinom.observed, multinom.null)
  all.deviance <- c(all.deviance, multinom.observed$deviance)
  all.deviance.factor <- c(all.deviance.factor, multinom.factor$deviance)
  all.df <- c(all.df, multinom.observed$edf)
  all.chi.squared <- c(all.chi.squared, linear.test$Chisq[2])
  all.p.values <- c(all.p.values, linear.test$`Pr(>Chisq)`[2])
  
  # CALCULATING PREDICTED VALUES
  
  coeff <- as.data.frame(coef(multinom.observed))
  predicted.0 <- matrix(1, nrow=tail(months.numeric,1)-months.numeric[1]+1,
                        ncol=length(rownames(coeff)), 
                        dimnames=list(c(months.numeric[1]:tail(months.numeric,1)), 
                                      rownames(coeff)))*c(months.numeric[1]:tail(months.numeric,1))
  predicted.1 <- t(apply(predicted.0, 1,
                         function(x) {exp(coeff$`(Intercept)` + coeff$MONTHS*x + coeff$MONTHS.SQUARED*x^2)}))

  # CALCULATING BASELINE
  
  div.by <- apply(predicted.1, 1, sum) + 1
  predicted.2 <- predicted.1/div.by
  baseline <- 1/div.by
  predicted.3 <- cbind(baseline, predicted.2)
  colnames(predicted.3)[1] = as.character(stats.categories[trait,])
  
  # CREATING PREDICTED DATA TABLE

  traits.predicted <- rbind(traits.predicted, t(predicted.3))
}

# CREATING DATA FRAMES

traits.predicted <- cbind(rownames(traits.predicted), traits.predicted)
colnames(traits.predicted)[1] <- "MODALITY"
traits.predicted <- join(modality, traits.predicted)
rownames(traits.predicted) <- traits.predicted$MODALITY
traits.predicted <- traits.predicted[-c(1,2)]
traits.predicted[is.na(traits.predicted)] <- 0
traits.predicted <- t(traits.predicted)
write.csv(as.matrix(traits.predicted), file=paste(path_out, 'FD_EPR_Predicted.csv', sep="/"))

# RANDOMIZATIONS

trait.rand.devs <- data.frame()
number.rands <- 1000

for (trait in rownames(stats.categories)){
  
  deviances.rand <- c()
  
  for (i in c(1:number.rands)){
    
    print(as.character(trait))
    print(i)
    
    # RANDOMIZING ABUNDANCE BY COLUMNS
    
    abund.rand <- abund[,sample(ncol(abund))]
    colnames(abund.rand) <- colnames(abund)
    print(abund.rand[1,])
    
    # CREATING DATA FRAME IN MULTINOM FORMAT
    
    epr.data.rand <- data.frame()
    for(month in months.character){
      multinom.format <- as.data.frame(lapply(traits, rep, abund.rand[,month]), check.names = FALSE)
      multinom.format <- cbind(rep(as.numeric(month), dim(multinom.format)[1]), multinom.format)
      epr.data.rand <- rbind(epr.data.rand, multinom.format)
    }
    
    colnames(epr.data.rand)[1] <- "MONTHS"
    epr.data.rand$MONTHS.SQUARED = as.numeric(epr.data.rand$MONTHS**2)
    
    print(dim(epr.data.rand))
    
    # CALCULATING DEVIANCE
    
    epr.data.rand$X <- droplevels(relevel(as.factor(epr.data.rand[,trait]), ref = as.character(stats.categories[trait,])))
    multinom.rand <- multinom(X~MONTHS+MONTHS.SQUARED, data = epr.data.rand, maxit=800)
    deviances.rand <- c(deviances.rand, multinom.rand$deviance)
  }
  trait.rand.devs <- rbind(trait.rand.devs, deviances.rand)
}

all.p.values <- rowSums(trait.rand.devs < all.deviance)/number.rands

deviance.table <- as.data.frame(as.numeric(format(round(all.deviance-all.deviance.factor, 0), nsmall = 0)), row.names = rownames(stats.categories))
colnames(deviance.table) <- "Deviance"

df.table <- as.data.frame(as.numeric(format(round(all.df, 0), nsmall = 0)), row.names = rownames(stats.categories))
colnames(df.table) <- "df"

chi.sq.table <- as.data.frame(as.numeric(format(round(all.chi.squared, 1), nsmall = 1)), row.names = rownames(stats.categories))
colnames(chi.sq.table) <- "Chisq"

p.val.table <- as.data.frame(as.numeric(format(round(all.p.values, 3), nsmall = 3)), row.names = rownames(stats.categories))
colnames(p.val.table) <- "p"

stats.categories <- cbind(stats.categories, deviance.table, df.table, chi.sq.table, p.val.table)

write.csv(stats.categories, file=paste(path_out, "FD_EPR_Statistics.csv", sep="/"), row.names = FALSE)
