#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 10:41:19 2017

@author: dawningliu
"""

import re
import os
import numpy as np
import matplotlib.pyplot as plt

# modules below are used in function "correctOnePressure", which was adapted from Paul E Ohno's spike removel script
import shutil
import pandas
from math import floor

# change to the directory of interest

# *********to run the script on your own computer, make sure to modify the path here************ 
os.chdir('/Users/dawningliu/Box Sync/research (YangdonglingLiu2020@u.northwestern.edu)/SFGdata/111317/isothermnew6_coc_1600')

# the function to generate the background spectrum by averaging data from multiple acquisitons
def aveBg(bgfname):   

    bgfhandle = open(bgfname,'r')
    bgfile = bgfhandle.readlines()
    del(bgfile[0],bgfile[1])
       
    replst = re.findall('[xX]([1-9]+)\.',bgfname)
    rep = replst[0]
    repnum = int(rep)
    
    wls_str = bgfile[0].split()[3:]
    wls = np.array([float(x) for x in wls_str])
    wns = 10**7/wls - 12500
    
    bglst = [[]]*repnum 
    for i in range(repnum):
        bglst[i] = [float(x) for x in bgfile[i+1].split()][2:]
    
    bgarray = np.array(bglst)
    
    bgaves = np.mean(bgarray,axis=0)
    
    np.savetxt('avebg.txt', bgaves)
    np.savetxt('wn.txt',wns)
    
    return 

# define a function to generate a background-subtracted sample spectrum (one DFG) at one pressure 
def aveOnePressure(fname,vp_voc):
    
    fhandle = open(fname,'r')
    file = fhandle.readlines()
    del(file[0],file[0],file[0])
    
    # to get the pressure from the filename
    
    # typical filename "He01_coc04_3x2.txt"
    # flow rate for Helium only  = 0.1 slpm, flow rate for Helium through cyclooctane = 0.4 slpm
    # 3 min acquisition time, 2 acquisitions
    # extract the two flow rates from the filename
    fratelst = re.findall('([0-9]+)_',fname)
    # flow rate for Helium only
    frateHe = int(fratelst[0])/10
    # flow rate for Helium through VOC (volatile organic compound)
    frateVoc = int(fratelst[1])/10
    # total flow rate
    frateTotal = frateHe + frateVoc
    
    # vp_voc is the SATURATED vapor pressure from user's input. See function isotherm line 194.
    # calculate the pressure of voc
    PVoc = (float(vp_voc))*(frateVoc/frateTotal)
    
    # extract the repetition number n from filename
    replst = re.findall('[xX](.+?)',fname)
    rep = replst[0]
    repnum = int(rep)
    
    
    lst = [[]]*repnum 
    for i in range(repnum):
        lst[i] = [float(x) for x in file[i].split()][2:]
    samplearray = np.array(lst)
    
    aves = np.mean(samplearray,axis=0)
    
    
    bgaves = np.loadtxt('avebg.txt')
       
    bgsubaves = aves - bgaves
    
    # save the background subtracted average as text file, indicating pressure in the filename
    np.savetxt('bgsubave_'+str(round(PVoc,2))+'.txt',bgsubaves)
    
    return

# This function is written by Paul E. Ohno, with my slight modifications for it to fit my project.
# it is to detect and correct non-meaningful spikes in the spectrum for each pressure
def correctOnePressure(fname):
    
    counts = np.loadtxt(fname)
    wns = np.loadtxt('wn.txt')
    
    # figure 1: plot the original spectrum
    plt.figure(1,figsize=(6,6))
    plt.plot(wns,counts)
    plt.xlabel('wavenumber',fontsize = 14)
    plt.ylabel('SFG Intensity',fontsize = 14)    
    
    windowSize = 7
    medians = pandas.Series(counts).rolling(window = windowSize,center=True).median()
    
    #replace beginning and end nan with the first/last computed value 
    numRep = floor(windowSize/2)
    for i in range(numRep):
        medians[i] = medians[numRep]
        medians[len(medians)-i-1] = medians[len(medians)-numRep-1]
    
    differences = counts-medians
    threshold = 100
    
    #empty array to hold zeros, or one if point is a spike
    spike = np.zeros(len(differences),)
    for i in range(len(differences)):
        if differences[i] > threshold or differences[i] < 0-threshold:
            spike[i] = 1
            print("Spike found at point index",i,'wavenumber',round(wns[i],2),'in',fname)
            
    #if a peak is found
    if np.sum(spike) > 0:
        
        origFile = open(fname,'r') # This line is essential!!! (but why???)
        newFile = open("temp" + fname,"w")
        
        countsCORR = counts.copy()
    
        for i in range(len(spike)):
            singleLine = origFile.readline()
            
            if spike[i] == 1:
                
                for j in range(5):
                    if (i-1-j) < 0:
                        left = [] 
                        break
                    else:
                        if spike[i-1-j] == 0:
                            left = counts[i-1-j] 
                            break
                
                for j in range(5):
                    if (i+j+1) >= len(spike):
                        right = [] 
                        break
                    else:
                        if spike[i+1+j] == 0:
                            right = counts[i+1+j] 
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

        newFile.close()
        shutil.copy2(fname,'raw'+fname)
        os.rename("temp"+fname,fname)
        
        # in figure 1: plot the corrected spectrum on top of the original one, to demonstrate the spike removal
        plt.plot(wns,countsCORR)
    
    else:
        print("No spikes found for",fname)
    return
    
# define a function to plot isotherms: sqrt. peak area vs pressure
def isotherm():
    
    # an input function will generate a pop-up question to let us enter vapor pressure for specific VOCs 
    # For cyclooctane (coc), vapor pressure = 4.6 torr
    vp_voc = input('Enter vapor pressure for this compound:')
    
    # go through the directory, detect background and sample data files, and work them up
    filenames = os.listdir()
    for filename in filenames:
        if re.search('.txt$',filename): # to make sure the filename ends with '.txt'
            if 'bg' in filename:
                aveBg(filename)
    for filename in filenames:
        if re.search('.txt$',filename): # to make sure the filename ends with '.txt'
            if 'He' in filename:
                aveOnePressure(filename,vp_voc)
    
    # correct spikes for each pressure, and overwrite the bgsubave files
    # original files became "raw..."
    morefilenames = os.listdir()
    for filename in morefilenames:
        if 'bgsubave' in filename:
            if 'raw' not in filename:
              correctOnePressure(filename)
            
    andmorefilenames = os.listdir()
    
    # create empty arrays for voc pressures and peak areas
    PVocs = np.array([])
    areasRaw = np.array([])
    
    for filename in andmorefilenames:
        if 'bgsubave' in filename:
            if not 'raw' in filename: # make sure it's after spike removal
                
                # fill the voc pressure array
                PVoclst = re.findall('_(.+)\.txt',filename) # extract the pressure from filename
                PVoc = float(PVoclst[0])
                PVocs = np.append(PVocs,PVoc)
                # fill the peak area array
                array = np.loadtxt(filename)
                area = np.sum(array)
                areasRaw = np.append(areasRaw,area)
    
    # subtract the smallest peak area from all the peak areas so that peak area = 0 when voc pressure = 0;
    # take square root of all peak areas
    sqrtareas = (areasRaw - np.amin(areasRaw))**0.5
    # normalize the array to (0,1)
    sqrtareasNorm = sqrtareas/np.amax(sqrtareas)

    # figure 2: plot normalized sqrt peak areas vs voc pressures
    plt.figure(2,figsize=(6,6))
    # dash lines connecting green round markers
    plt.plot(PVocs,sqrtareasNorm,'--go')
    plt.xlabel('Pressure',fontsize = 14)
    # subscript SFG in axis label
    plt.ylabel('$E_{SFG}$',fontsize = 14)
    
    # save both arrays for voc pressures and normalized sqrt peak areas as separate text files
    np.savetxt('pressures.txt',PVocs)
    np.savetxt('sqrtareasNorm.txt',sqrtareasNorm)

    return               

# run it!
isotherm()