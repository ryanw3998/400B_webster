#!/usr/bin/env python
# coding: utf-8

# In[11]:
'''Ryan Webster
ASTR 400B
Homework 2'''

import numpy as np
import astropy.units as u

def Read(filename): #defining function
    #Opening file
    file = open(filename,'r') #opening file
    
    #Extracting time
    line1 = file.readline() #This reads a line
    label, value = line1.split() #This deliniates the line from above to get time
    time = float(value)*u.Myr #assigning time a unit of mega year
    
    #Extracting number of particles
    line2 = file.readline() #same as above
    label, value = line2.split()
    num_part = float(value)
    
    file.close() #closing file
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #generating numpy array from remaining data
    
    return(time,num_part,data) #returning values of interest
    

