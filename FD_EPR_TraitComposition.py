#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FUNCTIONAL TRAIT COMPOSITION

Created on Mon Jun 11 14:23:09 2018

@author: laurendykman

Requires three input data tables:
    Abundance = species as rows, sites as columns. Contains a row of temperatures and a row of
        months associated with each sample. Temperature row must be named "Temperature". Months
        row must be named "Months". Names are not case-sensitive. Temperature and months must 
        be type int.
    Traits = species as rows, traits as columns. Cells contain modalities within a given trait.
        All modalities must be true modalities in data frame Modality. Case sensitive, space sensitive.
    Modality = table with traits in one column, corresponding modalities in another. Mdality names
        must exactly match the values in the Traits data frame. This script will run diagnostics
        to make sure all traits are present.
"""

# IMPORT FUNCTIONS AND PACKAGES

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# IMPORT FD_EPR MODULE 

import FD_EPR_TraitAnalysis as fd

# WORKING DIRECTORY

path = '/Users/laurendykman/Desktop/github/FD_EPR' # Enter your working directory
path_out = path + '/output'
os.chdir(path)

# IMPORT DATA

abund_raw = fd.Format("DataS1.csv")
trait = fd.Format("DataS2_BCO-DMO_Dykman.csv")
modality = fd.Format("FD_EPR_Modality.csv")
predicted = fd.Format(path_out + "/" + "FD_EPR_Predicted.csv")

# FORMATTING ABUNDANCE DATA

abund = abund_raw.drop(["Unnamed: 1"], axis=1)
abund = abund.drop(["Temperature"], axis=0)
index = abund.loc["Months Post-Eruption"] == "Pre-Eruption"
abund[abund.columns[index]] = abund[abund.columns[index]].apply(pd.to_numeric, errors='coerce', axis=1)
abund.replace(np.nan, -10, inplace=True)
abund = abund.groupby(abund.loc["Months Post-Eruption"], axis=1).sum()
months = list(abund.columns)
abund = abund.drop(["Months Post-Eruption"], axis=0)

# FORMATTING TRAIT DATA

index = trait.columns.str.contains('CITATION', case=False)
trait = trait.drop(trait.columns[index], axis=1) # Removing citation columns
trait = trait.drop("AphiaID", axis=1) # Removing Taxon Aphia ID column
trait = trait.drop("scientificNameID", axis=1)

# CHECK SPECIES LISTS FOR TRAIT AND ABUNDANCE TABLES ARE IDENTICAL

fd.CheckDataMatch(abund, trait)
fd.CheckModalities(trait, modality)

# PERFORMING MATRIX OPERATIONS TO GET RELATIVE ABUNDANCE OF EACH TRAIT MODALITY

trait_dict = fd.GetTraitDict(modality) # Using the modality table, generates a dictionary of modalities per trait.
trait_binary = fd.GetBinaryTraits(trait, modality) # Creates a binary table with "1" if an organism expresses a modality, and "0" otherwise.
trait_abund = fd.GetTraitAbundance(trait_binary, abund, modality) # Combines binary modality data with species abundance to calculate modality abundace.
trait_abund.to_csv(path_out + '/' + 'FD_EPR_TraitAbund.csv') # Uncomment the following line to save trait abundance as a csv file.

samples_to_plot = list(trait_abund.index) # Creating the list of all samples to plot for observed values.
samples_to_plot_fits = list(predicted.columns) # Creating a list of samples to plot for predicted values.

# CREATE A SEPARATE CWM PLOT FOR EACH TRAIT
# Change last command to True to save figures

for key in trait_dict:
    plt.figure(figsize=(10,5), dpi=100)
    fd.PlotTraitComposition(trait_abund, samples=None, modalities=trait_dict[key], fit=predicted, title=key.upper(), xlabel="Month", ylabel='Relative Abundance', filename = "FD_EPR_TRAIT_COMP", output_path = path_out, save = True, addline = True)
    
    
    
    