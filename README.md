# FD_EPR
A functional analysis of community-level change in a hydrothermal vent community recovering from volcanic disturbance.

This GitHub package contains data and scripts used to run analyses and produce figures for the publication: Dykman, L.N., Beaulieu, S.E., Solow, A.R., Mills, S.W., Mullineaux, L.S. (2021) Functional traits provide new insight into recovery and succession at deep-sea hydrothermal vents. Ecology. Accepted 2021-03-15.

Software Versions:

R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin17.0 (64-bit)

Python 3.8.1 (default, Jan  8 2020, 16:15:59) 
IPython 7.19.0 - An enhanced Interactive Python.

Data and Metadata Tables

DataS1.csv is related to Supplemental DataS1: [CountsPerSpeciesPerSandwich_EPR.csv] in Dykman et al. (2021) and is a subset of samples published as a dataset with additional metadata in BCO-DMO repository:

Mullineaux, L. (2020) Counts of colonists collected from colonization plates at the East Pacific Rise (EPR) deep-sea vents (1998-2017). Biological and Chemical Oceanography Data Management Office (BCO-DMO). (Version 2) Version Date 2020-08-31 doi:10.26008/1912/bco-dmo.733173.2

DataS2_BCO-DMO_Dykman.csv is related to Supplemental DataS1: [ModalitiesPerTraitPerSpecies_EPR.csv] in Dykman et al. (2021) and is published as a dataset with additional metadata in BCO-DMO repository:

Dykman, L.N. (2021) Functional traits of colonists collected from colonization surfaces at the East Pacific Rise (EPR) deep-sea vents from 1998-2017. Biological and Chemical Oceanography Data Management Office (BCO-DMO). (Version 1) Version Date 2020-03-21 doi:10.26008/1912/bco-dmo.844993.1

MetadataS1 is related to Supplemental DataS1: [ModalitiesPerTraitDefinitions_EPR.csv] in Dykman et al. (2021), and provides definitions and additional information for trait modalities used in the study. FD_EPR_Modality.csv is similar to this supplementary file, but is a reformatted version that is used in analyses. It is formatted to facilitate indexing and matrix operations for some of the scripts that run analyses, and includes guild names.

Workflow:

1)	Run FD_EPR_GowerDistance.R:

It is important to run this script first because it defines the functional guilds of species, which are used in all future analyses.

Inputs:
DataS2_BCO-DMO_Dykman.csv - a table of modalities per trait per species, with species as rows and modalities as columns.

Functions: Calculates Gower Distance. Can handle numerical, ordinal, or categorical data types. Changes data type of four columns from categorical to ordered categorical. Can handle missing values. Creates a relatedness dendrogram for 58 species based on similarities in modality assignments for 8 functional traits. Assigns functional guilds (designated as capital letters) to species. At the beginning of the script, removes the final column “FUNCTIONAL_GUILD” from DataS2_BCO-DMO_Dykman.csv to exclude guild assignments from dissimilarity calculations. At the end of the script, newly-calculated guild assignments are appended as the new final column. Overwrites DataS2_BCO-DMO_Dykman.csv with new functional guilds.

Outputs:
FD_EPR_DENDROGRAM_YYYY-MM-DD.pdf - a figure showing a dissimilarity dendrogram for all species based on dissimilarity in their modalities for 8 traits. This is Figure 3a in Dykman et al. (2021), and is created in an output subfolder.
Note: overwrites DataS2_BCO-DMO_Dykman.csv in main folder with updated final column “FUNCTIONAL_GUILD”.

2)	Run FD_EPR_Multinom.R

This script produces the statistics presented in the paper for changes in trait modality and functional guild composition over time. It uses multinomial logistic regression in the package “nnet” to calculate coefficients and deviance for the relative changes in modalities within each trait and for functional guilds as a function of time. This model is quadratic in time to allow for non-linear trends. P-values are calculated by 1,000 randomizations, which involves repeatedly randomizing columns while keeping observation times fixed. Due to the large number of randomizations, this script may take ~24hrs to run.

Inputs:
DataS1.csv
DataS2_BCO-DMO_Dykman.csv
FD_EPR_Modality.csv

Outputs:
FD_EPR_Predicted.csv - table of predicted values calculated by multinomial logistic regression. These are later plotted as fit lines in Figures 2 and 3b in Dykman et al (2021).
FD_EPR_Statistics.csv - an output file with deviance, degrees of freedom, chi-squared values, and p-values for each of the eight traits and for functional guilds. This is presented in Dykman et al. (2021) as Supplemental AppendixS5.

3)	Run FD_EPR_TraitAnalysis.py

FD_EPR_TraitAnalysis.py includes necessary packages and custom definitions, functions, and source code for all analyses conducted in Python. Run this before running any other Python script to assure these dependencies are loaded.

4)	Run FD_EPR_TraitComposition.py

This script imports FD_EPR_TraitAnalysis.py as fd. Figures are plotted using the function def PlotTraitComposition( ) in FD_EPR_TraitAnalysis.py.

Inputs:
DataS1.csv
DataS2_BCO-DMO_Dykman.csv
FD_EPR_Modality.csv
FD_EPR_Predicted.csv - from outputs folder, produced by FD_EPR_Multinom.R

Outputs: Figures 2 and 3b from Dykman et al. (2021)
FD_EPR_TRAIT_COMP_EXTERNAL PROTECTION_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_FEEDING METHOD_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_FUNCTIONAL GUILD_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_HABITAT COMPLEXITY_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_LARVAL DEVELOPMENT_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_MAXIMUM ADULT BODY SIZE_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_RELATIVE ADULT MOBILITY_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_REPRODUCTIVE TYPE_YYYY-MM-DD.pdf
FD_EPR_TRAIT_COMP_TROPHIC MODE_YYYY-MM-DD.pdf

Optional: Uncomment line 74 to save a trait abundance table, which has months since the 2006 eruption as rows and trait modalities as columns, as FD_EPR_TraitAbund.csv. Note: in this table, the “Pre-Eruption” sample is saved as timepoint “-10.” This is for plotting purposes only – the pre-eruption sample was taken in 1998.

5)	FD_EPR_Diversity.R

Uses the package FD to calculate Rao’s Quadratic Entropy (RaoQ) and Functional Richness (FRic). Run this R diversity file (FD_EPR_Diversity.R) before running the diversity file in Python (FD_EPR_Diversity.py), because the former uses outputs from the latter.

Inputs:
DataS1.csv
DataS2_BCO-DMO_Dykman.csv

Outputs:
FD_EPR_RaoQ.csv – a table of RaoQ values calculated by FD, which are plotted in python in a later step using the script FD_EPR_Diversity.py.
FD_EPR_FRic.csv – a table of FRic values calculated by FD, which are plotted in python in a later step using the script FD_EPR_Diversity.py.

6)	FD_EPR_Diversity.py

Calculates Hill diversity of order 1 for both species and guilds. Hill diversity is calculated by the function def Hill( ) in FD_EPR_TraitAnalysis.py. Imports Rao’s Quadratic Entropy (RaoQ) and Functional Richness (FRic) as tables output by FD_EPR_Diversity.R in Step 5. Statistics are calculated using the function def Randomize( ) and best fit lines are produced using the function def PolynomialFit( ).

Inputs:
DataS1.csv
DataS2_BCO-DMO_Dykman.csv
FD_EPR_Modality.csv
FD_EPR_RaoQ.csv - from outputs folder, written by FD_EPR_Diversity.R
FD_EPR_FRic.csv - from outputs folder, written by FD_EPR_Diversity.R

Outputs:
Figure 4 in Dykman et al. (2021) and AppendixS4
FD_EPR_DIVERSITY_SPECIES_YYYY-MM-DD.pdf
FD_EPR_DIVERSITY_GUILD_YYYY-MM-DD.pdf
FD_EPR_DIVERSITY_RAO_YYYY-MM-DD.pdf
FD_EPR_DIVERSITY_FRIC_YYYY-MM-DD.pdf
