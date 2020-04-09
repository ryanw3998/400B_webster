

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

# stream_array_x = [] 
# stream_array_y = [] 
# stream_array_z = [] 

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
        M31_disk_mass = ComponentMass(filenameM31,2)*1e12
        M33_disk_mass = ComponentMass(filenameM33,2)*1e12

        #Storing the COM pos and vel
        POS_MW = COM_MW.COM_P(0.1,5.0)
        VEL_MW = COM_MW.COM_V(POS_MW[0],POS_MW[1],POS_MW[2])

        #Storing the COM pos and vel
        POS_M31 = COM_M31.COM_P(0.1,5.0)
        VEL_M31 = COM_M31.COM_V(POS_M31[0],POS_M31[1],POS_M31[2])

        #Storing the COM pos and vel
        POS_M33 = COM_M33.COM_P(0.1,4.0)
        VEL_M33 = COM_M33.COM_V(POS_M33[0],POS_M33[1],POS_M33[2])

        #Centering M31
        X_M33_M31 = POS_M33[0] - POS_M31[0]
        Y_M33_M31 = POS_M33[1] - POS_M31[1]
        Z_M33_M31 = POS_M33[2] - POS_M31[2]

        #Finding the magnitude of the distance between M33 and M31
        R_M33_M31 = np.sqrt(X_M33_M31**2 + Y_M33_M31**2 + Z_M33_M31**2)

        #Calculating Jacobi Radius
        M33_r_jacobi =  R_M33_M31*(M33_disk_mass/(2*(M31_disk_mass)))**(1/3)
        print(M33_r_jacobi)

        #Centering M33 particles in array to apply 
        #Jacobi radius conditional

        # print(COM_M33.x[0],POS_M33[0])
        M33_part_x = COM_M33.x[:] - POS_M33[0].value
        M33_part_y = COM_M33.y[:] - POS_M33[1].value
        M33_part_z = COM_M33.z[:] - POS_M33[2].value
        particle_rad = np.sqrt(M33_part_x**2 + M33_part_y**2 + M33_part_z**2)
        print(len((particle_rad)))
        
        # determine the index for those particles within the Jacobi radius                                                
        index_jacobi = np.where(particle_rad >= M33_r_jacobi.value)

        # Applying condiational
        stream_x = COM_M33.x[:][index_jacobi]
        stream_y = COM_M33.y[:][index_jacobi]
        stream_z = COM_M33.z[:][index_jacobi]

        print(len(stream_x),len(stream_y),len(stream_z))

        # Plotting x and y star positions
        # plt.plot(stream_x,stream_y,'o')
        # plt.show()

        # stream_array_x.append(stream_x)
        # stream_array_y.append(stream_y)
        # stream_array_z.append(stream_z)

