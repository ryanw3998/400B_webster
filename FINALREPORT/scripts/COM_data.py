

"""
Ryan Webster
ASTR 400B 
Code that calculates COM positions and velocities
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

"""Script to calculate LG positional information throughout simulation.
This data is used in the velocity dispertion analysis. Whenever I need positional
or velocity data of one object wrt another, I would use this script. I uncomment the objects
I want to calculate and change the variable names accordingly. Messy but quick."""

start = 0
end = 800
n = 5
snap_ids = np.arange(start,end,n)

COM_pos_data = np.zeros((len(snap_ids),7))

#Looping over all low res simulation files
for  i, snap_id in enumerate(snap_ids):
    #composing the data filename 
    ilbl = '000' + str(snap_ids[i])
    ilbl = ilbl[-3:]
    print(ilbl)
    filenameMW = master_path+"%s_"%('MW')+"VLowRes/"+"%s_"%('MW') + ilbl + '.txt'
    filenameM31 = master_path+"%s_"%('M31')+"VLowRes/"+"%s_"%('M31') + ilbl + '.txt' 
    filenameM33 = master_path+"%s_"%('M33')+"VLowRes/"+"%s_"%('M33') + ilbl + '.txt' 
    
    # Initializing CenterOfMass classes for all galaxies
    # COM_MW = com2.CenterOfMass(filenameMW,2)
    # POS_MW = COM_MW.COM_P(0.1,5.0) # This POS array is COM of disk particles
    # VEL_MW = COM_MW.COM_V(POS_MW[0],POS_MW[1],POS_MW[2])

    # COM_M31 = com2.CenterOfMass(filenameM31,2)
    # POS_M31 = COM_M31.COM_P(0.1,5.0) # This POS array is COM of disk particles
    # VEL_M31 = COM_M31.COM_V(POS_M31[0],POS_M31[1],POS_M31[2])

    COM_M33 = com2.CenterOfMass(filenameM33,2)
    POS_M33 = COM_M33.COM_P(0.1,4.0) # This POS array is COM of disk particles
    VEL_M33 = COM_M33.COM_V(POS_M33[0],POS_M33[1],POS_M33[2])

    # # Creating arrays to find MW M31 COM
    # MW_M31_mass = np.concatenate((COM_MW.m0,COM_M31.m0))
    # MW_M31_pos_x = np.concatenate((COM_MW.x0,COM_M31.x0))
    # MW_M31_pos_y = np.concatenate((COM_MW.y0,COM_M31.y0))
    # MW_M31_pos_z = np.concatenate((COM_MW.z0,COM_M31.z0))
    # x_cm = (np.sum(MW_M31_mass*MW_M31_pos_x)/np.sum(MW_M31_mass))*u.kpc
    # y_cm = (np.sum(MW_M31_mass*MW_M31_pos_y)/np.sum(MW_M31_mass))*u.kpc
    # z_cm = (np.sum(MW_M31_mass*MW_M31_pos_z)/np.sum(MW_M31_mass))*u.kpc
    # POS_MWM31 = [x_cm,y_cm,z_cm]

    # Subtracting one object from another
    COM_pos_data[i,0] = ilbl
    COM_pos_data[i,1] = POS_M33[0].value# - POS_M31[0].value
    COM_pos_data[i,2] = POS_M33[1].value# - POS_M31[1].value
    COM_pos_data[i,3] = POS_M33[2].value# - POS_M31[2].value
    COM_pos_data[i,4] = VEL_M33[0].value# - VEL_M31[0].value
    COM_pos_data[i,5] = VEL_M33[1].value# - VEL_M31[1].value
    COM_pos_data[i,6] = VEL_M33[2].value# - VEL_M31[2].value

    # # This is for pos, vel info as it is the simulation (wrt MW)
    # COM_pos_data[i,1] = POS_M31[0].value
    # COM_pos_data[i,2] = POS_M31[1].value
    # COM_pos_data[i,3] = POS_M31[2].value
    # COM_pos_data[i,1] = VEL_M31[0].value
    # COM_pos_data[i,2] = VEL_M31[1].value
    # COM_pos_data[i,3] = VEL_M31[2].value

# print(COM_pos_data)
    
path_COM_data = '{}M33_data.txt'.format(master_path)
np.savetxt(path_COM_data, COM_pos_data, fmt = "%11.4f"*7, comments='#',
            header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                    .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))