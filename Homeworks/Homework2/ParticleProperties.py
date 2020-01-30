
'''Ryan Webster
ASTR 400B
Homework 2'''


import numpy as np
import astropy.units as u
from ReadFile import Read

filename = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/MW_000.txt' #defining file name

def ParticleInfo(filename,particle_type,particle_num): #defining function
    time, num_part, data = Read(filename) #Reading in file

    particle = data[data['type'] == particle_type][particle_num] #Finding particle based on inputs particle_type and particle_num

    mag_dist = float(np.around(np.sqrt((particle[2]**2)+(particle[3]**2)+(particle[4]**2)),3))*u.kpc #calculating magnitude of distance
    mag_vel = float(np.around(np.sqrt((particle[5]**2)+(particle[6]**2)+(particle[7]**2)),3))*(u.km/u.s) #calculating magnitude of velocity
    mass = float(particle[1])*10**10*u.M_sun #extracting particle mass and converting to M_sun

    return(mag_dist,mag_vel,mass) #returning values

#Testing functions
print('Homework 2 Question 5')
print('3D distance = ', ParticleInfo(filename,2,99)[0])#This returns the 3D distance, 3D velocity and mass for Quesiton 5 Parts 1-3 of HW2
print('3D velocity = ', ParticleInfo(filename,2,99)[1])
print('Mass = ', ParticleInfo(filename,2,99)[2])
distance_lyr = np.around(ParticleInfo(filename,2,99)[0].to(u.lyr),3)
print('3D distance (lyrs) = ',distance_lyr) #Converting 3D distance to lyrs for Questoin 5 Part 4 of HW2

