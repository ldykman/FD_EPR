#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FUNCTIONAL DIVERSITY

Created on Tue Jun 26 15:47:40 2018

@author: laurendykman
"""

# IMPORT FD_EPR MODULE

import FD_EPR_TraitAnalysis as fd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import os

# WORKING DIRECTORY

path = '/Users/laurendykman/Desktop/github/FD_EPR' # Enter your working directory
path_out = path + '/output'
os.chdir(path)

# IMPORTING DATA

abund_raw = fd.Format("DataS1.csv")
trait = fd.Format("DataS2_BCO-DMO_Dykman.csv")
modality = fd.Format("FD_EPR_Modality.csv")
raoQ = fd.Format(path_out + "/" + "FD_EPR_RaoQ.csv")
FRic = fd.Format(path_out + "/" + "FD_EPR_FRic.csv")

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
trait = trait.drop("scientificNameID", axis=1) # Removing scientific name ID column

# CHECK SPECIES LISTS FOR TRAIT AND ABUNDANCE TABLES ARE IDENTICAL

fd.CheckDataMatch(abund, trait)
fd.CheckModalities(trait, modality)

# CALCULATING TRAIT ABUNDANCE

trait_dict = fd.GetTraitDict(modality) # Using the modality table, generates a dictionary of modalities per trait.
trait_binary = fd.GetBinaryTraits(trait, modality) # Creates a binary table with "1" if an organism expresses a modality, and "0" otherwise.
trait_abund = fd.GetTraitAbundance(trait_binary, abund, modality) # Combines binary modality data with species abundance to calculate modality abundace.

# PLOTTING SPECIES DIVERSITY

hill_species = fd.Hill1(abund)

# Calculating quadratic fit and statistics

fitrange = list(range(1,int(months[-1])))
fit_species = fd.PolynomialFit(months[1:], hill_species[1:], 2)
predicted_species = np.polyval(fit_species['coeff'], fitrange)
pval_species = fd.Randomize(x = months[1:], y = hill_species[1:], rsq = fit_species['rsq'], deg = 2, n = 1000)

# Setting the plot

plt.figure(dpi=100)
plt.plot(fitrange, predicted_species, color= 'k')
a = plt.scatter(months[1:], hill_species[1:], color = 'k', marker = 'o', facecolors='none', 
               label = "rsq=" + str(round(float(fit_species['rsq']), 2)) + ", p=" + str(round(float(pval_species), 3)))
plt.axvline(x=0, color = 'k')
plt.scatter(-10, hill_species[0], color = 'k', marker = 'o')
plt.legend(handles = [a])
plt.text(-13, 11, "Pre")
plt.xlabel("Month", fontsize=15)
plt.ylabel("Hill Diversity", fontsize=15)
plt.title("Species Diversity", fontsize=15)
plt.savefig(path_out + "/" + "FD_EPR_DIVERSITY_SPECIES_"+str(datetime.date.today())+".pdf")

# PLOTTING GUILD DIVERSITY

guild = trait_abund[trait_dict['FUNCTIONAL GUILD']].transpose() # Subsets, only uses functional group column.
hill_guild = fd.Hill1(guild)

# Calculating quadratic fit and statistics

fit_guild = fd.PolynomialFit(months[1:], hill_guild[1:], 2)
predicted_guild = np.polyval(fit_guild['coeff'], fitrange)
pval_guild = fd.Randomize(months[1:], hill_guild[1:], fit_guild['rsq'], 2, 1000)

# Setting the plot

plt.figure(dpi=100)
plt.plot(fitrange, predicted_guild, color= 'k')
b = plt.scatter(months[1:], hill_guild[1:], color = 'k', marker = 'o', facecolors='none', 
               label="rsq=" + str(round(fit_guild['rsq'], 2)) + ", p=" + str(round(float(pval_guild), 3)))
plt.axvline(x=0, color = 'k')
plt.scatter(-10, hill_guild[0], color = 'k', marker = 'o')
plt.legend(handles = [b])
plt.text(-13, 3.7, "Pre")
plt.xlabel("Month", fontsize=15)
plt.ylabel("Hill Diversity", fontsize=15)
plt.title("Guild Diversity", fontsize=15)
plt.savefig(path_out + "/" + "FD_EPR_DIVERSITY_GUILD_"+str(datetime.date.today())+".pdf")

# CALCULATING INDICES OF FUNCTIONAL DIVERSITY

rao = np.array(raoQ["x"])
FRic = np.array(FRic["x"])

# PLOTTING RAOQ

plt.figure(dpi=100)
fitrange2 = list(range(int(months[2]),int(months[-1])))
fit_rao = fd.PolynomialFit(months[1:], rao[1:], 2)
fit2_rao = fd.PolynomialFit(months[2:], rao[2:], 2)
predicted_rao = np.polyval(fit_rao['coeff'], fitrange)
predicted2_rao = np.polyval(fit2_rao['coeff'], fitrange2)
pval_rao = fd.Randomize(months[1:], rao[1:], fit_rao['rsq'], 2, 1000)
pval2_rao = fd.Randomize(months[2:], rao[2:], fit2_rao['rsq'], 2, 1000)
plt.plot(fitrange, predicted_rao, color= 'k')
plt.plot(fitrange2, predicted2_rao, color= 'k', linestyle=':')
c = plt.scatter(months[1:], rao[1:], color = 'k', marker = 'o', facecolors='none',
               label="RaoQ (rsq=" + str(round(fit_rao['rsq'], 2)) + ", p=" + str(round(float(pval_rao), 2)) + ")")
plt.axvline(x=0, color = 'k')
plt.scatter(-10, rao[0], color = 'k', marker = 'o')
plt.legend(labels = ("rsq=" + str(round(fit_rao['rsq'], 2)) + ", p=" + str(round(float(pval_rao), 3)), "rsq=" + str(round(fit2_rao['rsq'], 2)) + ", p=" + str(round(float(pval2_rao), 3))))
plt.text(-13, 0.023, "Pre")
plt.xlabel("Month", fontsize=15)
plt.ylabel("RaoQ", fontsize=15)
plt.title("Rao's Quadratic Entropy", fontsize=15)
plt.savefig(path_out + "/" + "FD_EPR_DIVERSITY_RAO_"+str(datetime.date.today())+".pdf")

# PLOTTING FRIC

plt.figure(dpi=100)
fit_FRic = fd.PolynomialFit(months[1:], FRic[1:], 2)
predicted_FRic = np.polyval(fit_FRic['coeff'], fitrange)
pval_FRic = fd.Randomize(months[1:], FRic[1:], fit_FRic['rsq'], 2, 1000)
plt.plot(fitrange, predicted_FRic, color= 'k')
f = plt.scatter(months[1:], FRic[1:], color = 'k', marker = 'o', facecolors='none',
               label="rsq=" + str(round(fit_FRic['rsq'], 2)) + ", p=" + str(round(float(pval_FRic), 2)))
plt.axvline(x=0, color = 'k')
plt.scatter(-10, FRic[0], color = 'k', marker = 'o')
plt.legend(handles=[f])
plt.text(-13, 0.04, "Pre")
plt.xlabel("Month", fontsize=15)
plt.ylabel("FRic", fontsize=15)
plt.title("Functional Richness", fontsize=15)
plt.savefig(path_out + "/" + "FD_EPR_DIVERSITY_FRIC_"+str(datetime.date.today())+".pdf")








