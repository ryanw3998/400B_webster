

"""
Ryan Webster
ASTR 400B 
Code that finds M33 stream particles
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
master_path = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/research project/'

"""This script is used to find M33 stream particles over the 
course of the MW and M31 merger. It calculates 3 different Jacobi radii of M33 with
at every snap shot in the simulation. The three radii are between: M31 and M33, 
MW and M33, and the MW/M31 COM and M33. It uses the lowest of these three radii
to select M33's stellar stream stars and saves a .txt file containing this information."""

start = 0
end = 100
n = 5
snap_ids = np.arange(start,end,n)

# Bootstrapping initial Jacobi radius with aprroximate M31/M33 distance at
# the beginning of the simultaion
r_jacobi = 210*u.kpc
MW_M33_rj = r_jacobi
M31_M33_rj = r_jacobi
MWM33COM_M33_rj = r_jacobi

# Calculating the number of data points to expect for r_jacobi_array
npoins = np.rint((end-start)/n).astype(int)

# This is an array that will contain the computed jacobi radii
# It is used to plot the minimum jacobi radius at each point
# in the simulation.
r_jacobi_array = np.zeros((npoins,4))

#Looping over all low res simulation files
for  i, snap_id in enumerate(snap_ids):
    #composing the data filename 
    ilbl = '000' + str(snap_ids[i])
    ilbl = ilbl[-3:]
    print(ilbl)
    filenameMW = master_path+"%s_"%('MW')+"VLowRes/"+"%s_"%('MW') + ilbl + '.txt'
    filenameM31 = master_path+"%s_"%('M31')+"VLowRes/"+"%s_"%('M31') + ilbl + '.txt' 
    filenameM33 = master_path+"%s_"%('M33')+"VLowRes/"+"%s_"%('M33') + ilbl + '.txt' 
    
    #Initializing CenterOfMass classes for all galaxies
    COM_MW = com2.CenterOfMass(filenameMW,2)
    POS_MW = COM_MW.COM_P(0.1,5.0) # This POS array is COM of disk particles

    COM_M31 = com2.CenterOfMass(filenameM31,2)
    POS_M31 = COM_M31.COM_P(0.1,5.0) # This POS array is COM of disk particles

    COM_M33 = com2.CenterOfMass(filenameM33,2)
    POS_M33 = COM_M33.COM_P(0.1,4.0) # This POS array is COM of disk particles

    # Calculating Jacobi Radius for MW and M33    
    MW_M33_sep = com2.separation(POS_MW,POS_M33)
    MW_mencl = COM_MW.MassEnclosedTotal([MW_M33_sep.value], 5.0)
    MW_M33_mencl = COM_M33.MassEnclosedTotal([r_jacobi.value], 4.0)

    MW_M33_rj = com2.JacobiRadius(MW_M33_sep,MW_M33_mencl,MW_mencl)
    print('MW_M33_rj',MW_M33_rj)

    # Calculating Jacobi Radius for M31 and M33
    M31_M33_sep = com2.separation(POS_M31,POS_M33)
    M31_mencl = COM_M31.MassEnclosedTotal([M31_M33_sep.value], 5.0)
    M31_M33_mencl = COM_M33.MassEnclosedTotal([r_jacobi.value], 4.0)

    M31_M33_rj = com2.JacobiRadius(M31_M33_sep,M31_M33_mencl,M31_mencl)
    print('M31_M33_rj',M31_M33_rj)

    #Calculating Jacobi Raidus for M31/MW center of mass and M33
    # First, I am calculating the COM of MW and M31
    MW_M31_mass = np.concatenate((COM_MW.m0,COM_M31.m0))
    MW_M31_pos_x = np.concatenate((COM_MW.x0,COM_M31.x0))
    MW_M31_pos_y = np.concatenate((COM_MW.y0,COM_M31.y0))
    MW_M31_pos_z = np.concatenate((COM_MW.z0,COM_M31.z0))
    x_cm = (np.sum(MW_M31_mass*MW_M31_pos_x)/np.sum(MW_M31_mass))*u.kpc
    y_cm = (np.sum(MW_M31_mass*MW_M31_pos_y)/np.sum(MW_M31_mass))*u.kpc
    z_cm = (np.sum(MW_M31_mass*MW_M31_pos_z)/np.sum(MW_M31_mass))*u.kpc
    POS_MWM31 = [x_cm,y_cm,z_cm]

    # Calculating the distance between MW/M31 COM and M33
    MWM33COM_M33_sep = com2.separation(POS_MWM31,POS_M33)

    # Now, I need to calculate the mass enclosed of MW and M31 to the distance of MWM33COM_M33_sep
    # First, I find the distance from the MW to M33
    MW_wrt_COM_x = COM_MW.x0 - x_cm.value #Note .x0 is x position for all paritcles, same for .y0 and .z0
    MW_wrt_COM_y = COM_MW.y0 - y_cm.value
    MW_wrt_COM_z = COM_MW.z0 - z_cm.value
    MW_wrt_COM_r = np.sqrt(MW_wrt_COM_x**2 + MW_wrt_COM_y**2 + MW_wrt_COM_z**2)
    # Now, I locate MW particles which are within MWM33COM_M33_sep and sum their masses
    MWencl = np.sum(COM_MW.m0[np.where(MW_wrt_COM_r<MWM33COM_M33_sep.value)])

    # Same thing as above, but for M31 now
    M31_wrt_COM_x = COM_M31.x0 - x_cm.value
    M31_wrt_COM_y = COM_M31.y0 - y_cm.value
    M31_wrt_COM_z = COM_M31.z0 - z_cm.value
    M31_wrt_COM_r = np.sqrt(M31_wrt_COM_x**2 + M31_wrt_COM_y**2 + M31_wrt_COM_z**2)
    M31encl = np.sum(COM_M31.m0[np.where(M31_wrt_COM_r<MWM33COM_M33_sep.value)])

    # Adding together the enclosed masses
    M_encl_tot = MWencl + M31encl

    # Calculating Jacobi radius of M33 with MW/M31 COM
    MWM33COM_M33_mencl = COM_M33.MassEnclosedTotal([r_jacobi.value], 4.0)
    MWM33COM_M33_rj = com2.JacobiRadius(MWM33COM_M33_sep,MWM33COM_M33_mencl,M_encl_tot*1e10*u.Msun)
    print('MWM33COM_M33_rj',MWM33COM_M33_rj)

    # Finding smallest jacobi radius to select stream stars
    r_jacobi = np.min((MW_M33_rj,M31_M33_rj,MWM33COM_M33_rj))*u.kpc
    print('Minimum Jacobi Radius', r_jacobi)

    # Saving jacobi radii to r_jacobi_array for later analysis
    r_jacobi_array[i,0] = MW_M33_rj[0].value
    r_jacobi_array[i,1] = M31_M33_rj[0].value
    r_jacobi_array[i,2] = MWM33COM_M33_rj[0].value
    r_jacobi_array[i,3] = r_jacobi.value

    # Centering M33 particles in array to apply Jacobi radius conditional     
    M33_part_x = COM_M33.x - POS_M33[0].value # disk particles
    M33_part_y = COM_M33.y - POS_M33[1].value
    M33_part_z = COM_M33.z - POS_M33[2].value
    particle_rad = np.sqrt(M33_part_x**2 + M33_part_y**2 + M33_part_z**2)

    # determine the index for those particles within the Jacobi radius                                                
    index_jacobi = np.where(particle_rad > r_jacobi.value)
    # Applying condiational
    stream_x = COM_M33.x[index_jacobi]
    stream_y = COM_M33.y[index_jacobi]
    stream_z = COM_M33.z[index_jacobi]
    stream_vx = COM_M33.vx[index_jacobi]
    stream_vy = COM_M33.vy[index_jacobi]
    stream_vz = COM_M33.vz[index_jacobi]
    stream_m = COM_M33.m[index_jacobi]
    snap = np.full(len(index_jacobi[0]),snap_id)

    # At the beginning of the simulation there are no stream particles
    # This conditional will save a txt file of zeros in place of an empty streams
    # array
    if len(stream_x) == 0:
        data = np.stack(([snap_id],[0],[0],[0],[0],[0],[0],[0]),axis=1)
    # If there are stream particles, then a text file is saved containing the
    # snap id, and particle mass, position and velocity.
    else:
        data = np.stack((snap,stream_m,stream_x,stream_y,stream_z,\
            stream_vx,stream_vy,stream_vz),axis=1)

    # Saving stream particles at each snapshot
    print("Number of Stream Particles = {}".format(len(index_jacobi[0])))
    savename = 'stream_{}.txt'.format(ilbl)
    path = '{}stream_data/{}'.format(master_path,savename)
    np.savetxt(path, data, fmt = "%11.4f"*8, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

    # Centering simulation on Andromeda
    xstream_cen = stream_x-POS_M31[0].value
    ystream_cen = stream_y-POS_M31[1].value
    xM31_cen = POS_M31[0].value - POS_M31[0].value
    yM31_cen = POS_M31[1].value - POS_M31[1].value
    xM33_cen = POS_M33[0].value - POS_M31[0].value
    yM33_cen = POS_M33[1].value - POS_M31[1].value
    xMW_cen = POS_MW[0].value - POS_M31[0].value
    yMW_cen = POS_MW[1].value - POS_M31[1].value

    #Creating live view of streams
    plt.close()
    plt.figure(figsize=(7,6))
    plt.gca().set_aspect(1)
    plt.xlabel('x (kpc)')
    plt.ylabel('y (kpc)')
    plt.xlim(-400,400)
    plt.ylim(-400,400)
    plt.title('M33 Streams (Snap ID: {})'.format(snap_id))
    plt.plot(xstream_cen,ystream_cen,'go',markersize=.1)
    plt.plot(xM33_cen,yM33_cen,'ro',markersize=1)
    plt.plot(xM31_cen,yM31_cen,'bo',markersize=1)
    plt.plot(xMW_cen,yMW_cen,'yo',markersize=1)
    plt.savefig('{}stream_video_better/{}.png'.format(master_path,ilbl))
    plt.pause(0.001)

#Saving jacobi radii
# np.savetxt('{}jacobi_radii_min.txt'.format(master_path),r_jacobi_array)

# Plotting all jacobi radii as a function of time
start = 0
end = 800
n = 5
r_jacobi_array = np.loadtxt('{}jacobi_radii_min.txt'.format(master_path))
label = ['MW','M31','MW/M31 COM','Radius Used']
plt.title('M33 Jacobi Radius for Local Group')
plt.ylabel('Jacobi Radius (kpc)')
plt.xlabel('Snap Shot')
plt.plot(np.arange(start,end,n),r_jacobi_array)
plt.legend(label)
plt.savefig('{}plots/M33_jacobi_radii_min.png'.format(master_path))
plt.show()