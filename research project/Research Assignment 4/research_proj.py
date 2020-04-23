

"""
Ryan Webster
ASTR 400B 
Draft of code for research project
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from IPython.display import Latex
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass
import MassProfile as mp

"""This is my first attempt at my research code.
I am planning on studying the stellar streams of M33 
by finding particles outside the Jacobi Radius of the 
galaxy over the course of the simulation. Once 
I find those particles, I can move on to other parts of
my project. For example, studying their velocity gradients
and velocity dispersions. 

I would like to plot the particles outside the Jacobi 
radius. 

The idea for the method invloving the Jacobi radius 
is from Dr. Besla's comments to my research proposal.
She asked how I was planning on finding stellar streams
around M33, and suggested using the Jacobi radius."""

""" This code is similar to code used in OrbitCOM.py
It will loop through snap shots of the merger simulation
and create a CenterOfMass instance for each galaxy. 
I will use COMP to find the distance of M33 from M31.
The jacobi radius is given by R(M_sat/(2*M_host))^(1/3).
I plug in the M31/M33 separation for R, and use ComponetMass()
to find M_sat and M_host. Once I have the Jacobi radius, I will
find stars in M33's disk that are outside of said radius.
For this assignment, I will simply plot those points."""


start = 0
end = 800
n = 5
snap_ids = np.arange(start,end,n)
# streams = np.zeros([8,1])

# print(np.shape(streams))
r_jacobi = 210*u.kpc

for  i, snap_id in enumerate(snap_ids):# loop over files
    #composing the data filename 
    ilbl = '000' + str(snap_ids[i])
    ilbl = ilbl[-3:]
    filenameMW = "%s_"%('MW')+"VLowRes/"+"%s_"%('MW') + ilbl + '.txt'
    filenameM31 = "%s_"%('M31')+"VLowRes/"+"%s_"%('M31') + ilbl + '.txt' 
    filenameM33 = "%s_"%('M33')+"VLowRes/"+"%s_"%('M33') + ilbl + '.txt' 
    
    #Initializing CenterOfMass class
    COM_MW = CenterOfMass(filenameMW,2)# Uses disk particles
    COM_M31 = CenterOfMass(filenameM31,2)# Uses disk particles
    COM_M33 = CenterOfMass(filenameM33,2)# Uses disk particles  
      

    #Storing the COM pos and vel MW
    POS_MW = COM_MW.COM_P(0.1,5.0)
    VEL_MW = COM_MW.COM_V(POS_MW[0],POS_MW[1],POS_MW[2])     
    #Storing the COM pos and vel M31
    POS_M31 = COM_M31.COM_P(0.1,5.0)
    VEL_M31 = COM_M31.COM_V(POS_M31[0],POS_M31[1],POS_M31[2])          
    #Storing the COM pos and vel M33
    POS_M33 = COM_M33.COM_P(0.1,4.0)
    VEL_M33 = COM_M33.COM_V(POS_M33[0],POS_M33[1],POS_M33[2])     

    # Centering M31
    X_M33_M31 = POS_M33[0] - POS_M31[0]
    Y_M33_M31 = POS_M33[1] - POS_M31[1]
    Z_M33_M31 = POS_M33[2] - POS_M31[2]        

    #Finding the magnitude of the distance between M33 and M31
    R_M33_M31 = np.sqrt(X_M33_M31**2 + Y_M33_M31**2 + Z_M33_M31**2)
    M31_mencl = COM_M31.MassEnclosedTotal([R_M33_M31.value], 5.0)
    M33_mencl = COM_M33.MassEnclosedTotal([r_jacobi.value], 4.0)

      
    #Calculating Jacobi Radius
    r_jacobi = R_M33_M31*(M33_mencl/(2*(M31_mencl)))**(1/3) #assumes isothermal sphere
    print("M31 M33 distance = {}".format(R_M33_M31))
    print("Jacobi Radius = {}".format(r_jacobi[0]))

    # #Centering M33 particles in array to apply Jacobi radius conditional     
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

    if len(stream_x) == 0:
        data = np.stack(([snap_id],[0],[0],[0],[0],[0],[0],[0]),axis=1)
    else:
        data = np.stack((snap,stream_m,stream_x,stream_y,stream_z,\
            stream_vx,stream_vy,stream_vz),axis=1)

    print("Number of Stream Particles = {}".format(len(index_jacobi[0])))
    savename = 'stream_{}.txt'.format(ilbl)
    path = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/research project/Code/streams/{}'.format(savename)
    np.savetxt(path, data, fmt = "%11.4f"*8, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
  
    # # Plotting x and y star positions
    print(POS_M31[0].value)
    print(stream_x)
    xstream_cen = stream_x-POS_M31[0].value
    ystream_cen = stream_y-POS_M31[1].value
    xM31_cen = POS_M31[0].value - POS_M31[0].value
    yM31_cen = POS_M31[1].value - POS_M31[1].value
    xMW_cen = POS_MW[0].value - POS_M31[0].value
    yMW_cen = POS_MW[1].value - POS_M31[1].value

    #Creating live view of streams
    plt.clf()
    plt.xlabel('kpc')
    plt.ylabel('kpc')
    plt.xlim(-500,500)
    plt.ylim(-500,500)
    plt.title('M33 Streams (Snap ID: {})'.format(snap_id))
    plt.gca().set_aspect(1)
    plt.plot(xstream_cen,ystream_cen,'go',markersize=.1)
    plt.plot(X_M33_M31,Y_M33_M31,'ro',markersize=1)
    plt.plot(xM31_cen,yM31_cen,'bo',markersize=1)
    plt.plot(xMW_cen,yMW_cen,'yo',markersize=1)
    plt.savefig('streams_video2/{}.png'.format(ilbl))
    # plt.pause(0.001)

    # plotting velocity dispersion ???
    plt.plot(stream_x,stream_vx,'o',ms=1)
