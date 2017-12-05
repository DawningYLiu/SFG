#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:41:13 2017

@author: dawningliu
"""

import re
import os
import numpy as np
import matplotlib.pyplot as plt

# modules below are used in function "correctFinalSpec", which was adapted from Paul E Ohno's spike removel script
# import shutil for moving and copying files
import shutil
# import pandas for rolling median function
import pandas
from math import floor

# change to the directory of interest
# *********to run the script on your own computer, make sure to modify the path here************ 
os.chdir('/Users/dawningliu/Box Sync/research (YangdonglingLiu2020@u.northwestern.edu)/SFGdata/100817/PS4T3_092914_9-12AM_S2')

# Define a function to generate the background spectrum by averaging data from multiple acquisitons
def aveBg(bgfname):   
    
    #import the background file
    bgfhandle = open(bgfname,'r')
    bgfile = bgfhandle.readlines()
    
    # delete the irrelevant lines - the first line: index numbers
    # delete the irrelevant lines - the third line: strings "Frame" "Strip".
    # "bgfile[1]" operates on the file with "bgfile[0]" executed, i.e. first line deleted.
    del(bgfile[0],bgfile[1]) 
    
    
    # typical filename "1300_4x3.txt" means DFG=1300, 4 min per acquisition, 3 acquisitions.
    # Using regular expression to extract the repetition number n of data acquisitoin.
    # n is the same for both bg and sample
    replst = re.findall('[xX]([1-9]+)\.',bgfname) 
    rep = replst[0]
    repnum = int(rep)
    # repnum should be equal to (len(bgfile)+1)
    
    
    # extract the wavelength list from the background data file
    # exclude the first three non-meaningful entries: '"', '"', '"wavelegnth"'
    wls_str = bgfile[0].split()[3:]
    
    # create numpy arrays instead of sticking with lists!
    wls = np.array([float(x) for x in wls_str])
    
    # convert the wavelength (nm) of SFG light to the wavenumber (1/cm) of IR light
    wns = 10**7/wls - 12500
    
    
    # make a list of n empty lists
    bglst = [[]]*repnum 
   
    # extract each line of counts as a list. n lines in total.
    for i in range(repnum):
        bglst[i] = [float(x) for x in bgfile[i+1].split()][2:] # exclude the first two non-meaningful entries
    
    # make an ndarray for the list of lists.
    bgarray = np.array(bglst)
    
    # calculate the mean of the n background acquisitions 
    bgaves = np.mean(bgarray,axis=0)
    
    # save wavenumbers and averaged background intensity as text files
    np.savetxt('avebg.txt', bgaves)
    np.savetxt('wn.txt',wns)
    
    return 

# define a function to generate a background-subtracted sample spectrum at one DFG position
def aveOneDFG(fname):
    
    # import sample files
    fhandle = open(fname,'r')
    file = fhandle.readlines()
    
    # delete the first three irrelevant lines.
    # the second line is the wavelength list, which I already got from the background file
    del(file[0],file[0],file[0])
    
    
    # typical filename "1300_4x3.txt" means DFG=1300, 4 min per acquisition, 3 acquisitions.
    # extract the DFG position from the filename
    dfglst = re.findall('([0-9]+)_',fname)
    dfg = dfglst[0]
    # we don't need DFG number for calculation, so it can remain as a string
    
    # extract the repetition number n, same for both bg and sample
    replst = re.findall('[xX](.+?)',fname)
    rep = replst[0]
    repnum = int(rep)
    
    
    # make a list of n empty lists
    lst = [[]]*repnum 
    
    # extract each line of counts as a list. n lines in total.
    for i in range(repnum):
        lst[i] = [float(x) for x in file[i].split()][2:]
    
    # make an ndarray for the list of lists.
    samplearray = np.array(lst)
    
    # calculate the mean of the n sample data acquisitions 
    aves = np.mean(samplearray,axis=0)
    
    
    # load background data from the previously generated file
    bgaves = np.loadtxt('avebg.txt')
       
    # subtract the background mean from the sample mean
    bgsubaves = aves - bgaves
    
    
    # save the background subtracted average as text file, indicating dfg position in the filename
    np.savetxt('bgsubave'+dfg+'.txt',bgsubaves)
    
    return


# define a function to generate a full spectrum by summing up spectra at all DFGs
def sumDFGs():
    
    # go through the directory, detect background and sample data files, and work them up
    filenames = os.listdir()
    for filename in filenames:
        if re.search('.txt$',filename): # to make sure the filename ends with '.txt'
            if 'bg' in filename:
                aveBg(filename)
    for filename in filenames:
        if re.search('.txt$',filename): # to make sure the filename ends with '.txt'
            if not 'bg' in filename:
                aveOneDFG(filename)
    
    
    wns = np.loadtxt('wn.txt')
    length = len(wns)
    
    
    # after generating worked-up data files for each DFG, sum them up to generate a full spectrum
    morefilenames = os.listdir()
    sumcounts = np.zeros(length)
    for filename in morefilenames:
        if 'bgsubave' in filename:
            bgsubaves = np.loadtxt(filename)
            sumcounts += bgsubaves
    
    # save the full spectrum data as text file
    np.savetxt('fullsfg.txt',sumcounts)
    
    return

# This function is written by Paul E. Ohno, with my slight modifications for it to fit my project.
# define a function to detect and correct non-meaningful spikes in the final full spectrum
# This function also plots the corrected spectrum
def correctFinalSpec(fname): 
    
    # load the full spectrum and the wavenumbers
    counts = np.loadtxt('fullsfg.txt')
    wns = np.loadtxt('wn.txt')
    
    # figure 1: plot the original spectrum
    plt.figure(1,figsize=(6,6))
    plt.plot(wns,counts,'b')
    # as a convention in IR and SFG/SHG spectroscopy, the x axis for wavenumbers is inverted
    plt.gca().invert_xaxis()
    plt.xlabel('wavenumber',fontsize = 14)
    plt.ylabel('SFG Intensity',fontsize = 14)
    
    
    # scan through the whole spectrum and find medians for every windowSize points
    windowSize = 7
    medians = pandas.Series(counts).rolling(window = windowSize,center=True).median()
    
    #number of nan at beginning and end to replace
    numRep = floor(windowSize/2)
    #replace beginning and end nan with the first/last computed value
    for i in range(numRep):
        medians[i] = medians[numRep]
        medians[len(medians)-i-1] = medians[len(medians)-numRep-1]
    
    #find difference of each point with the median of its window
    differences = counts-medians
    
    #threshold, past which if it is further from median it will sense that it is a spike
    threshold = 100
    
    #empty array to hold zeros, or one if point is a spike
    spike = np.zeros(len(differences),)
    for i in range(len(differences)):
        # negative spikes could come from spikes in background data
        if differences[i] > threshold or differences[i] < 0-threshold:
            spike[i] = 1
            print("Spike found at point index",i,'wavenumber',round(wns[i],2))
            
    #if a peak is found
    if np.sum(spike) > 0:
        
        #read in datafile
        origFile = open(fname,'r') 
        
        #create new file to put modified
        newFile = open("temp" + fname,"w")
        
        #create copy for new corrected array
        countsCORR = counts.copy()
    
        for i in range(len(spike)):
            singleLine = origFile.readline()
            
            #if the point needs to be replaced
            if spike[i] == 1:
                
                #check up to five points to the left for the edge or for an ok point
                for j in range(5):
                    if (i-1-j) < 0:
                        left = [] #if its edge only take from right point
                        break
                    else:
                        if spike[i-1-j] == 0:
                            left = counts[i-1-j] #or get the first acceptable point
                            break
                
                #check up to five points to the right for the edge or for an ok point        
                for j in range(5):
                    if (i+j+1) >= len(spike):
                        right = [] #if its edge only take from the left point
                        break
                    else:
                        if spike[i+1+j] == 0:
                            right = counts[i+1+j] #or get the first acceptable point
                            break
                
                #get the average of the two or the value if its only one
                tempValArray = np.array([])
                tempValArray = np.append(tempValArray,left)
                tempValArray = np.append(tempValArray,right)
                ave = tempValArray.mean()
                countsCORR[i] = ave
                
                #get line from original file, modify
                singleLine = str(ave) + '\n'
           
            #write original or modified line
            newFile.write(singleLine)

        # close new file
        newFile.close()
        
        # rename original file to preSpikeCorr 
        shutil.copy2(fname,"preSpikeCorr.txt")
        
        # rename the new temp file to the original name
        # the corrected full spectrum data in this file can be copy-pasted to Igor easily
        os.rename("temp"+fname,fname)
        
        # in figure 1: plot the corrected spectrum on top of the original one, to demonstrate the spike removal
        plt.plot(wns,countsCORR,'g')
        
        # figure 2: plot only corrected value, in a standardized way 
        plt.figure(2,figsize=(6,6))
        plt.plot(wns,countsCORR,'g')
        plt.gca().invert_xaxis()
        plt.xlabel('wavenumber',fontsize = 14)
        plt.ylabel('SFG Intensity',fontsize = 14)
    
    else:
        print("No spikes found")
    return    

# define a function that combines all the functions, to generate a full spectrum from multiple DFGs
def fullSpecWorkup():
    sumDFGs()
    correctFinalSpec('fullsfg.txt')


# run it!
fullSpecWorkup()

            
    
    
    


        
    




