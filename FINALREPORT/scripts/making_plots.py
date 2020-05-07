"""
Ryan Webster
ASTR 400B 
Code to plot data
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
import glob
from analyze_stream_particles import generate_video,velocity_dispersion
master_path = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/research project/'


"""Script to plot data found in analyze_stream_particles. Its a bit messy."""

# Importing stream particle data, selecting only snap 465 for velocity dispersion analysis
stream_pos_465 = np.sort(glob.glob("{}stream_data/stream_465.txt".format(master_path)))
M33_data = np.genfromtxt('{}M33_wrt_M31_data.txt'.format(master_path),dtype=None,names=True)

# Defining simulation start, end and interval
start = 0
end = 800
n = 5

# Plotting all jacobi radii as a function of time
r_jacobi_array = np.loadtxt('{}jacobi_radii_min.txt'.format(master_path))
label = ['MW','M31','MW/M31 COM','Radius Used']
plt.title('M33 Jacobi Radius for Local Group')
plt.ylabel('Jacobi Radius (kpc)')
plt.xlabel('Snap Shot')
plt.plot(np.arange(start,end,n),r_jacobi_array)
plt.legend(label)
plt.savefig('{}plots/M33_jacobi_radii_min.png'.format(master_path))
plt.show()

# Generating movie of M33 Streams during LG merger
generate_video(start,end,n)

# # Generating plot of velocity dispersion
# velocity_dispersion(stream_pos_465,0,150,n_shells=10,leading=True)

# # Centering stream particles on M33
# stream_xcen = stream_pos_465['x'] - M33_data['x'][93]
# stream_ycen = stream_pos_465['y'] - M33_data['y'][93]
# stream_zcen = stream_pos_465['z'] - M33_data['z'][93]

# # Indexing stream particles to isolate leading tail
# index1 = ((stream_ycen>=0) & (stream_ycen<=150))
# index2 = (stream_xcen<=25)
# stream_xcen_sel = stream_xcen[index1 & index2]
# stream_ycen_sel = stream_ycen[index1 & index2]

# # Loading in velocity dispersion data
# disp_data = np.genfromtxt('{}dispersion_data.txt'.format(master_path),dtype=None,names=True)

# # Creating subplots
# fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(10,4),gridspec_kw={'width_ratios': [1.5,1],'height_ratios': [1]})    

# # Plotting stream velocity dispersion
# ax[0].set_title("Stream Velocity Dispersion Snap: 465")
# ax[0].set_xlabel('Radial Distance (kpc)')
# ax[0].set_ylabel("Velocity Dispersion (km/s)")
# ax[0].plot(disp_data['radii'],disp_data['sigma_r'],linewidth=2,color='red',label=r'$\sigma_r$')
# ax[0].plot(disp_data['radii'],disp_data['sigma_t'],linewidth=2,color='blue',label=r'$\sigma_t$')
# ax[0].legend(prop={'size':12})

# # Plotting M33 streams at snap 465, highlighting streams used in dispersion analysis
# ax[1].set_title('M33 Streams Snap: 465')
# ax[1].set_xlabel('x (kpc)')
# ax[1].set_ylabel('y (kpc)')
# ax[1].set_aspect('equal', 'box')
# ax[1].set_xlim(-250,250)
# ax[1].set_ylim(-250,250)
# ax[1].plot(stream_xcen,stream_ycen,'go',markersize=.5,label='Stream Particles')
# ax[1].plot(stream_xcen_sel,stream_ycen_sel,'o',color='orange',markersize=.5,label='Leading Particles')
# ax[1].plot(0,0,'ko',markersize=1,label='M33')
# ax[1].legend(markerscale=5,prop={'size':9})

# # Saving plot
# plt.tight_layout()
# plt.savefig("{}plots/dispersion_and_snap.png".format(master_path),dpi=300)
# plt.show()