#Homework 6
#Ryan Webster
#Collaborators: Sean Cunningham, Madison Walder, Jimmy Lilly

#import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import matplotlib.pyplot as plt
import matplotlib
from ReadFile import Read
from CenterOfMass2 import CenterOfMass



def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
    galaxy: Name of desired galaxy
    start: Starting Snapshot
    end: Ending SnapShot
    n: Interval when to return COM
          
    returns:
    File of galactic position and velocity over start and end times 
    """
    
    #composing filename for output
    fileout = "Orbit_{}.txt".format(galaxy)
    
    #setting tolerance and VolDec
    if galaxy == 'M33':
        VolDec = 4.0
        delta = 0.1
    else:
        VolDec = 5.0
        delta = 0.1
    #print(galaxy,VolDec,delta)
    #generating the snapshot id sequence 
    snap_ids = np.arange(start,end,n)
    #initializing the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])
    
    for  i, snap_id in enumerate(snap_ids):# loop over files
        #composing the data filename 
        ilbl = '000' + str(snap_ids[i])
        ilbl = ilbl[-3:]
        filename = "%s_"%(galaxy)+"VLowRes/"+"%s_"%(galaxy) + ilbl + '.txt' #I kept my VLowRes files in a separate folder, 
        #to keep my Homework6 directory from being polluted with .txt files

        #Initializing CenterOfMass class
        COM = CenterOfMass(filename,2)# Uses disk particles
        #Storing the COM pos and vel
        POS = COM.COM_P(delta,VolDec)
        VEL = COM.COM_V(POS[0],POS[1],POS[2])
    
        #storting t, x, y, z, vx, vy, vz in obrit array
        orbit[i]= COM.time.value/1000, *tuple(POS.value), *tuple(VEL.value)
        
        #print(snap_id)
        
    #Writing out data
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

#Running funtion
#OrbitCOM('MW', 0, 800, 5)
#OrbitCOM('M31', 0, 800, 5)
#OrbitCOM('M33', 0, 800, 5)

#Reading in data
MW_data = np.genfromtxt("Orbit_MW.txt",dtype=None,names=True)
M31_data = np.genfromtxt("Orbit_M31.txt",dtype=None,names=True)
M33_data = np.genfromtxt("Orbit_M33.txt",dtype=None,names=True)

#Function used to calculate magnitude of differences between vectors of different galaxies
def mag_diff(data1,data2):

    xdiff = data1['x']-data2['x'] #subtracting x positions
    ydiff = data1['y']-data2['y'] #"" y ""
    zdiff = data1['z']-data2['z'] #"" z ""

    vxdiff = data1['vx']-data2['vx'] #subtracting x velocities
    vydiff = data1['vy']-data2['vy'] #"" y ""
    vzdiff = data1['vz']-data2['vz'] #"" z ""

    r = np.sqrt(xdiff**2.0+ydiff**2.0+zdiff**2.0) #Calculating magnitude of position vectors
    vr = np.sqrt(vxdiff**2.0+vydiff**2.0+vzdiff**2.0) #"" velocity ""

    return(r,vr)

#Using mag_diff for MW and M31
MW_M31_r,MW_M31_vr = mag_diff(MW_data,M31_data)

#Using mag_diff for M33 and M31
M33_M31_r,M33_M31_vr = mag_diff(M33_data,M31_data)


# Plot the Orbit of the galaxies 
#################################

plt.figure(1)
plt.title("Local Group Separation (0-11 Gyr)")
plt.xlabel('Time (Gyr)')
plt.ylabel('Separation (kpc)')
plt.plot(MW_data['t'],MW_M31_r,label='MW M31 Separation')
plt.plot(MW_data['t'],M33_M31_r,label='M33 M31 Separation')
plt.semilogy()
plt.legend()
plt.show()
#plt.savefig('local_group_sep.png',dpi=350)


# Plot the orbital velocities of the galaxies 
#################################

plt.figure(2)
plt.title("Local Group Velocity (0-11 Gyr)")
plt.xlabel('Time (Gyr)')
plt.ylabel('Velocity (km/s)')
plt.plot(MW_data['t'],MW_M31_vr,label='M33 M31 Velocity')
plt.plot(MW_data['t'],M33_M31_vr,label='M33 M31 Velocity')
plt.semilogy()
plt.legend()
plt.show()
#plt.savefig('local_group_vel.png',dpi=350)


#********** Section 4 **********

#Question 1:

#Based on the graph of galaxy separation vs time, MW and M31 will have 3 close encounters.
#This is based off the fact that the graph gets close to 0 kpc 3 times over the course of
#the simulation.

#Question 2:

#It appears that when the separation between MW and M31 is a minimum, the velocity is a maximum.
#This is also seen for M33 and M31. 

#Question 3:

#MW and M31 appear to merge after 6.6 Gyr. After the merger, M33 appears to orbit MW and M31.