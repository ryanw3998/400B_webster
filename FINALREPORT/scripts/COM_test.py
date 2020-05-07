

"""
Ryan Webster
ASTR 400B 
Center of Mass test script
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from IPython.display import Latex
import CenterOfMass2 as com2
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass
import MassProfile as mp

"""This is script used for testing the COM_P function in 
CenterofMass2. I was interested in finding if using all galaxy particles
significantly changed the COM positons of the galaxy. It does not, 
so I just stuck with using disk particles to find COM."""

start = 0
end = 800
n = 5
snap_ids = np.arange(start,end,n)

# Bootstrapping initial Jacobi radius with aprroximate M31/M33 distance at
# the beginning of the simultaion
r_jacobi = 210*u.kpc
MW_M33_rj = r_jacobi
M31_M33_rj = r_jacobi
MWM33COM_M33_rj = r_jacobi

npoins = np.rint((end-start)/5).astype(int)
r_jacobi_array = np.zeros((npoins,4))

#Initializing arrays for plotting differences in galaxy COM calculations
# as a funcitn of snap shot.
POS_MW_array = []
POS_MW_test_array = []

POS_M31_array = []
POS_M31_test_array = []

POS_M33_array = []
POS_M33_test_array = []

for  i, snap_id in enumerate(snap_ids):# loop over files
    #composing the data filename 
    print(i)
    ilbl = '000' + str(snap_ids[i])
    ilbl = ilbl[-3:]
    filenameMW = "%s_"%('MW')+"VLowRes/"+"%s_"%('MW') + ilbl + '.txt'
    filenameM31 = "%s_"%('M31')+"VLowRes/"+"%s_"%('M31') + ilbl + '.txt' 
    filenameM33 = "%s_"%('M33')+"VLowRes/"+"%s_"%('M33') + ilbl + '.txt' 
    
    #Initializing CenterOfMass class
    COM_MW = com2.CenterOfMass(filenameMW,2)
    POS_MW = COM_MW.COM_P(0.1,5.0) # This POS array is COM of disk particles
    POS_MW_test = COM_MW.COM_P_test(0.1,5.0)# This POS array is COM of all particles
    POS_MW_array.append(POS_MW.value)
    POS_MW_test_array.append(POS_MW_test.value)

    COM_M31 = com2.CenterOfMass(filenameM31,2)
    POS_M31 = COM_M31.COM_P(0.1,5.0) # This POS array is COM of disk particles
    POS_M31_test = COM_M31.COM_P_test(0.1,5.0)# This POS array is COM of all particles
    POS_M31_array.append(POS_M31.value)
    POS_M31_test_array.append(POS_M31_test.value)

    COM_M33 = com2.CenterOfMass(filenameM33,2)
    POS_M33 = COM_M33.COM_P(0.1,4.0) # This POS array is COM of disk particles
    POS_M33_test = COM_M33.COM_P_test(0.1,4.0)# This POS array is COM of all particles
    POS_M33_array.append(POS_M33.value)
    POS_M33_test_array.append(POS_M33_test.value) 

#Converting to numpy arrays, they are much easier to index fo graphing.
POS_MW_array = np.array(POS_MW_array)
POS_MW_test_array = np.array(POS_MW_test_array)

POS_M31_array = np.array(POS_M31_array)
POS_M31_test_array = np.array(POS_M31_test_array)

POS_M33_array = np.array(POS_M33_array)
POS_M33_test_array = np.array(POS_M33_test_array)

#Ploting MW test
fig,ax = plt.subplots(3,figsize=(8,6), sharex= True)
fig.suptitle('MW: All Particles vs Disk Particles COM')
plt.xlabel('Snap Shot')
ax[0].set(ylabel='x (kpc)')
ax[1].set(ylabel='y (kpc)')
ax[2].set(ylabel='z (kpc)')
snaps_array = np.arange(start,end,n)
ax[0].plot(snaps_array,POS_MW_array[:,0])
ax[0].plot(snaps_array,POS_MW_test_array[:,0])
ax[1].plot(snaps_array,POS_MW_array[:,1])
ax[1].plot(snaps_array,POS_MW_test_array[:,1])
ax[2].plot(snaps_array,POS_MW_array[:,2])
ax[2].plot(snaps_array,POS_MW_test_array[:,2])
ax[0].tick_params(axis='both',direction='in',length=6)
ax[1].tick_params(axis='both',direction='in',length=6)
ax[2].tick_params(axis='both',direction='in',length=6)
plt.legend(('Disk Particles','All Particles'),bbox_to_anchor=(0., 3.45, 1., .102), loc='lower center',
           ncol=2, borderaxespad=0.)
plt.savefig('plots/COM_test_MW.png')
plt.show()

#Ploting M31 test
fig,ax = plt.subplots(3,figsize=(8,6), sharex= True)
fig.suptitle('M31: All Particles vs Disk Particles COM')
plt.xlabel('Snap Shot')
ax[0].set(ylabel='x (kpc)')
ax[1].set(ylabel='y (kpc)')
ax[2].set(ylabel='z (kpc)')
snaps_array = np.arange(start,end,n)
ax[0].plot(snaps_array,POS_M31_array[:,0])
ax[0].plot(snaps_array,POS_M31_test_array[:,0])
ax[1].plot(snaps_array,POS_M31_array[:,1])
ax[1].plot(snaps_array,POS_M31_test_array[:,1])
ax[2].plot(snaps_array,POS_M31_array[:,2])
ax[2].plot(snaps_array,POS_M31_test_array[:,2])
ax[0].tick_params(axis='both',direction='in',length=6)
ax[1].tick_params(axis='both',direction='in',length=6)
ax[2].tick_params(axis='both',direction='in',length=6)
plt.legend(('Disk Particles','All Particles'),bbox_to_anchor=(0., 3.45, 1., .102), loc='lower center',
           ncol=2, borderaxespad=0.)
plt.savefig('plots/COM_test_M31.png')
plt.show()

#Ploting M33 test
fig,ax = plt.subplots(3,figsize=(8,6), sharex= True)
fig.suptitle('M33: All Particles vs Disk Particles COM')
plt.xlabel('Snap Shot')
ax[0].set(ylabel='x (kpc)')
ax[1].set(ylabel='y (kpc)')
ax[2].set(ylabel='z (kpc)')
snaps_array = np.arange(start,end,n)
ax[0].plot(snaps_array,POS_M33_array[:,0])
ax[0].plot(snaps_array,POS_M33_test_array[:,0])
ax[1].plot(snaps_array,POS_M33_array[:,1])
ax[1].plot(snaps_array,POS_M33_test_array[:,1])
ax[2].plot(snaps_array,POS_M33_array[:,2])
ax[2].plot(snaps_array,POS_M33_test_array[:,2])
ax[0].tick_params(axis='both',direction='in',length=6)
ax[1].tick_params(axis='both',direction='in',length=6)
ax[2].tick_params(axis='both',direction='in',length=6)
plt.legend(('Disk Particles','All Particles'),bbox_to_anchor=(0., 3.45, 1., .102), loc='lower center',
           ncol=2, borderaxespad=0.)
plt.savefig('plots/COM_test_M33.png')
plt.show()