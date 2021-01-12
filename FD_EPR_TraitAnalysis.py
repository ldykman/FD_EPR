#pl.co#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 12:17:20 2018

@author: laurendykman

DEFINITIONS, FUNCTIONS, PACKAGES, AND SOURCE CODE FOR FUNCTIONAL TRAIT ANALYSIS

"""

### IMPORT LIBRARIES ###

import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

### FUNCTIONS FOR IMPORTING CHECKING AND MANIPULATING DATA SETS ###
   
def Format(filename):
    
    """
    Parameters
    ----------
    filename: input the name of your csv data table as a string.
    Data should be formatted with samples as rows and species as columns.

    Returns
    -------
    A pandas data frame in the correct format for analysis.
    """
    
    file = pd.read_csv(filename)
    file.columns = [i for i in file.columns]
    file.index = file.iloc[:,0]
    file = file.drop(file.columns[0], axis=1)
    file.index = [i for i in file.index]
    return(file)
    
def CheckDataMatch(abundance, traits):
    
    """
    Checks whether species in the abundance and traits data frames are identical.
    If not, prints a warning and lists the species in question. If all species match,
    prints "All species match." To pass this test, the species names in the abundance
    and traits data frames must be in the same order and spelled identically. Punctuation
    and capitalization count!
    
    Parameters
    ----------
    abundance: a data frame with abundance count data, with species as rows 
    and samples as columns.
    traits: a data frame with trait modlity information for each species, with
    species as rows and traits as columns.
    """
    
    names1 = list(traits.index)
    names2 = list(abundance.index)
    if names1 != names2:
        for name in names1:
            if name not in names2:
                print("Warning, the following species are not in abundance matrix:")
                print(name)
        for name in names2:
            if name not in names1:
                print("Warning, the following species are not in traits matrix:")
                print(name)
    else:
        print("All species match.")
        
def CheckModalities(traits, modalities):
    
    """
    This function checks all the modality assignments in the traits data frame for each
    species to assure they are valid modalities. If unaccepted modalities are located,
    this funciton will print the incorrect modality and indicate which species
    has the error. Spelling and punctuation count. Will also alert the user to NA values, 
    implying modalities were unassigned due to insufficient information. These are 
    accepted in analysis and NA values are ignored.
    
    Parameters
    ----------
    traits: a data frame with trait modlity information for each species, with
    species as rows and traits as columns.
    modalities: a data frame with modality options for every trait. This data frame
    has two columns, the first lists the trait name, the second lists the modality 
    names.
    """
    
    counter = 0
    for j in range(0, len(traits.columns)):
        for i in range(0,len(traits.index)):
            counter += 1
            if (traits.iloc[i,j] in list(modalities.iloc[:,0])) == False:
                print(str(traits.iloc[i,j])+" is not a real modality. This error is in "
                  +str(traits.index[i])+", "+str(traits.columns[j]))
    print(str(counter)+" cells scanned in the trait matrix.")

def GetTraitDict(modalities):
    
    """
    Parameters
    ----------
    modalities: a data frame with modality options for every trait. This data frame
    has two columns, the first lists the trait name, the second lists the modality 
    names.

    Returns
    -------
    trait_dict: a python dictionary with traits as keys and modalities as values
    generated from a data frame that contains traits in the first column and possible
    modalities in the second column.
    """
    
    trait_dict = {}
    for name in modalities.index:
        trait_dict[name] = list(modalities.loc[name, modalities.columns[0]])
    return trait_dict

def GetBinaryTraits(traits, modalities):
    
    """
    Parameters
    ----------
    traits: a data frame with trait modlity information for each species, with
    species as rows and traits as columns.
    modalities: a data frame with modality options for every trait. This data frame
    has two columns, the first lists the trait name, the second lists the modality 
    names.

    Returns
    -------
    binary_trait_data: a pandas data frame with species as rows and all possible 
    modalities as columns indicating whether a species expresses that modality, where
    a value of 1 indicates the species expresses the modality and a value of 0 indicates
    the species does not express the modality.
    """
    
    modalities_list = list(modalities.iloc[:,0])
    species_list = list(traits.index)
    species_dict = {}
    for species in species_list:
        species_modalities = [str(i) for i in traits.loc[species]]
        binary_list = []
        for modality in modalities_list:
            if str(modality) in species_modalities:
                binary_list.append(1)
            else:
                binary_list.append(0)
        species_dict[species] = binary_list
    trait_data_binary = pd.DataFrame.from_dict(species_dict, orient='index')
    trait_data_binary.columns = modalities_list
    return trait_data_binary

def GetTraitAbundance(trait_data_binary, abundance, modalities):
    
    """
    Parameters
    ----------
    trait_data_binary: a data frame of binary values indicating whether a specices (rows)
    expresses a modality (columns). This is a data frame generated from the function 
    GetBinaryTraits().
    abundance: a data frame with abundance count data, with species as rows 
    and samples as columns.
    modalities: a data frame with modality options for every trait. This data frame
    has two columns, the first lists the trait name, the second lists the modality 
    names.

    Returns
    -------
    trait_data_abundance: a data frame showing counts for all organisms expressing 
    a modality (columns) for each sample (rows). These counts are calculated by summing
    the counts of all species that express a given modality in a given sample.
    """
    
    modality_abundance_dict = {}
    for sample, abund_list in abundance.iteritems():
        modality_sum_list = []
        for modality, binary_list in trait_data_binary.iteritems():
            modality_sum = sum(np.array(abund_list)*np.array(binary_list))
            modality_sum_list.append(modality_sum)
        modality_abundance_dict[sample] = modality_sum_list
    trait_data_abundance = pd.DataFrame.from_dict(modality_abundance_dict, orient='index')
    trait_data_abundance.columns = trait_data_binary.columns
    return trait_data_abundance
    
def PlotTraitComposition(trait_data_abundance, samples = None, modalities = None, fit = None, title = '', xlabel = None, ylabel = None, filename = '', output_path = None, save = False, addline = False):
    
    """
    Parameters
    ----------
    trait_data_abundance: a data frame showing counts for all organisms expressing 
        a modality (columns) for each sample (rows). These counts are calculated by summing
        the counts of all species that express a given modality in a given sample. Generated
        by the function GetTraitAbundance().
    samples: a list of sample names that you wish to include in the figure. If no 
        sample names are provided, the default is None, and all samples in trait_data_abundance 
        will be plotted.
    modalities: a list of modalities that you wish to plot. Normally, in calling this function,
        you would input all possible modalities for a single trait at a time. If no modalities
        are provided, the default is None, and all modalities in trait_data_abundance will be
        plotted. This is not ideal, as plotting the relative abundance of modalitites for multiple
        traits at once is meaningless.
    fit: an optional data frame showing statistical fits which will be printed as lines on top
        of the bar plot. This data frame should have time points as rows and modalities as columns,
        the same as trait_data_abundance. The default is None, in which case the observed 
        trait composition data will be plotted without statistical fits.
    title: a title for your plot, input as a string. The default is ''.
    xlabel: a label for the x axis of your plot. The default is None.
    ylabel: a label for the y axis of your plot. The default is None.
    filename: a filename for your figure, input as a string. The default is ''.
    save: True or False, would you like to save your figure as a pdf to the current working directory.
        The default is False.
    addline: would you like to plot the vertical line that separates pre-eruption samples
        from post-eruption samples. The default is False.
    """
    
    # Subsetting Data
    
    if samples == None:
        samples = list(trait_data_abundance.index)
    if modalities == None:
        modalities = list(trait_data_abundance.columns)
    trait_data_subset = trait_data_abundance.loc[samples, modalities]
    
    # Making Data Proportional

    trait_data_proportion = np.divide(trait_data_subset.transpose(), trait_data_subset.sum(axis=1)).transpose()
    
    # Setting the Plot
    
    spaces_bar = trait_data_proportion.index # Getting spacing along the x axis from months since the 2006 eruption
    x_offset_bar = np.array([0,-1,2,0,0,0,0,0]) # Two of the early samples are so close in time the bars overlap, so they are shifted slightly to left and right on figure for visibility.
    y_offset_bar = np.zeros(len(spaces_bar)) # Setting a counter to plot bars as cumulative relative abundance on the y axis.
    width = 7 # Width of bars
    index = 0 # Sets a counter
    
    # Plotting Composition and Statistical Fit Data
    
    for modality in modalities:
        colors = plt.cm.coolwarm(np.linspace(0.2, 1, len(modalities))) # Creating a list of colors on the purple scale.
        # For bars of observed data
        plt.bar(spaces_bar + x_offset_bar, trait_data_proportion[modality], width, bottom=y_offset_bar, color=colors[index], label = modalities[index])
        y_offset_bar += trait_data_proportion[modality] # Adding to the cumulative abundance so bars will be stacked.
        index += 1 # Incrementing the counter by 1
    
    # Plotting statistical fits

    if fit is None:
        pass
    else:
        index = 0
        fit_subset = fit[modalities]
        fit_proportion = np.divide(fit_subset.transpose(), fit_subset.sum(axis=1)).transpose()
        spaces_fit = fit_proportion.index # Getting spacing along the x axis for the predicted fits
        y_offset_fit = np.zeros(len(spaces_fit))
        for modality in modalities:
            colors = plt.cm.coolwarm(np.linspace(0.2, 1, len(modalities))) # Creating a list of colors on the purple scale.
            y_offset_fit += np.array(fit_proportion[modality])
            plt.plot(spaces_fit, y_offset_fit, color=colors[index])
            index += 1 # Incrementing the counter by 1
    
    # The following script formats the number of columns in the legend so it fits well
    
    colnum = int(len(modalities))
    if colnum == 4 or 7 < colnum < 9: 
        colnum = 2
    if 4 < colnum < 8 or colnum >= 9:
        colnum = 3
        
    # Plot aesthetics specifications

    plt.legend(bbox_to_anchor=(0.5, -0.25), loc=9, ncol=colnum, borderaxespad=0.0, fontsize=15)
    plt.title(title, fontsize=30)
    plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    plt.text(-16,-0.07,"Pre",fontsize=20) # Indicates the pre-eruption sample
    plt.xticks([0,20,40,60,80,100,120,140], fontsize=20) # x ticks every 20 months
    plt.yticks(fontsize=20)
    plt.ylim(0,1)
    #plt.tight_layout()
    if addline == True: # This vertical line visually separates pre-eruption samples from post-eruption samples
        plt.axvline(x=0.5, color='k')
    
    # Saving Plot
    
    if save == True:
        plt.savefig(output_path + '/' + filename+'_'+title.upper()+'_'+str(datetime.date.today())+'.pdf', bbox_inches='tight')
    
# FUNCTIONAL DISTANCE BETWEEN SPECIES

def Difference(trait):
    diff_matrix = []
    for species1 in trait.index:
        traits1 = np.array(list(trait.loc[species1]))
        row = []
        for species2 in trait.index:
            traits2 = np.array(list(trait.loc[species2]))
            diff = len(traits1) - sum(traits1 == traits2)
            row.append(diff)
        diff_matrix.append(row)
    diff_matrix = pd.DataFrame(diff_matrix)
    diff_matrix.index = trait.index
    diff_matrix.columns = trait.index
    return(diff_matrix)

# DIVERSITY INDICES

def Hill1(abundance):
    p = np.divide(abundance, abundance.sum(axis=0))
    log_p = np.log(p)
    log_p[log_p == -np.inf] = 0
    shan = -(p * log_p).sum()
    hill = np.array(np.exp(shan))
    return hill

# STATISTICS

def PolynomialFit(x, y, degree):
        
        results = {}
        coeffs = np.polyfit(x, y, degree)
        # Polynomial Coefficients
        results['coeff'] = coeffs.tolist()
        # r-squared
        p = np.poly1d(coeffs)
        # fit values, and mean
        yhat = p(x)
        ybar = np.sum(y)/len(y)
        ssreg = np.sum((yhat-ybar)**2)
        sstot = np.sum((y - ybar)**2)
        results['rsq'] = ssreg / sstot
        return results
    
def Randomize(x, y, rsq, deg, n):
    x = np.array(x)
    y = np.array(y)      
    rsq_list_rand = np.array([])
    for i in range(0,n):
        #print(i)
        np.random.shuffle(y)
        #print(y)
        fit = PolynomialFit(x, y, 2)
        rsq_list_rand = np.append(rsq_list_rand, fit['rsq'])
    pval = np.array(rsq_list_rand >= rsq).sum()/len(rsq_list_rand)
    return pval