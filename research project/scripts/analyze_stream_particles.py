"""
Ryan Webster
ASTR 400B 
Code to analyze stream particles
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
master_path = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/research project/'


"""Script to analyze stream particles found in find_stream_particles.py"""

# Importing stream particle data
stream_files = np.sort(glob.glob("{}stream_data/stream_*.txt".format(master_path)))

# Importing positional data of all local group objects (positions are wrt M31)
M33_data = np.genfromtxt('{}M33_pos_data.txt'.format(master_path),dtype=None,names=True)
MW_data = np.genfromtxt('{}MW_pos_data.txt'.format(master_path),dtype=None,names=True)
M31_data = np.genfromtxt('{}M31_pos_data.txt'.format(master_path),dtype=None,names=True)
MMWM31COM_pos_data = np.genfromtxt('{}MWM31COM_pos_data.txt'.format(master_path),dtype=None,names=True)

# Opening stream data for all snap shots
lengths = []
all_snaps = []
for i,file in enumerate(stream_files):
    open_file = np.genfromtxt("{}".format(file),dtype=None,names=True)
    all_snaps.append(open_file)
    lengths.append(len(open_file.data))
all_snaps = np.array(all_snaps)
lengths = np.array(lengths)

def generate_video(start,end,n,plotlim=.3):
    """Function that generates a movie of M33 streams in xy, xz, and yz planes.
    The total orbit of M33 is represented as a solid black line, while M33 is 
    a black dot on that line. M31 is represented as a blue dot at the center
    of the plot. The MW is represented as a red dot. The M31 and Mw COM is 
    represented as a purple x. The main purpose of this plot is to see whether
    M33's streams trace its orbit.

    Note: I converted distances to Mpc to make smaller axis labels
    
    Inputs:
    start: starting snap of simulation
    end: ending snap of simulation
    n: spacing of snap shots
    plotlim: Size of plot (Mpc)
    
    Returns:
    Movie of M33 streams in xy, xz, and yz planes
    """

    # Arranging snap shots to be used
    snaps = np.arange(start,end,n)

    # Using all M33 positional data
    # This is used to visualize the entirety of M33's orbit
    M33_all_posx = (M33_data['x']-M31_data['x'])/1000
    M33_all_posy = (M33_data['y']-M31_data['y'])/1000
    M33_all_posz = (M33_data['z']-M31_data['z'])/1000
    #Looping through snap shots
    for i,snap_id in enumerate(snaps):
        ilbl = '000' + str(snaps[i])
        ilbl = ilbl[-3:]

        #Selecting stream data at specific snapshot
        streamx = (all_snaps[i]['x']-M31_data[i]['x'])/1000
        streamy = (all_snaps[i]['y']-M31_data[i]['y'])/1000
        streamz = (all_snaps[i]['z']-M31_data[i]['z'])/1000

        #Finding M33's position at a snap shot
        M33_curr_posx = (M33_data[i]['x']-M31_data[i]['x'])/1000
        M33_curr_posy = (M33_data[i]['y']-M31_data[i]['y'])/1000
        M33_curr_posz = (M33_data[i]['z']-M31_data[i]['z'])/1000

        #Finding the MW's position at a snap shot
        MW_posx = MW_data[i]['x']/1000
        MW_posy = MW_data[i]['y']/1000
        MW_posz = MW_data[i]['z']/1000

        #Finding M33/MW COM position at a snap shot
        MWM31_com_posx = MMWM31COM_pos_data[i]['x']/1000
        MWM31_com_posy = MMWM31COM_pos_data[i]['y']/1000
        MWM31_com_posz = MMWM31COM_pos_data[i]['z']/1000

        #Creating live view of streams
        plt.close()
        fig, ax = plt.subplots(nrows = 1,ncols= 3,figsize=(14,5.5))
        fig.suptitle('M33 Streams (Snap ID: {})'.format(snap_id),fontsize=15)

        #Setting axis limits
        ax[0].set(xlabel='x (Mpc)',ylabel='y (Mpc)',xlim=(-plotlim,plotlim),ylim=(-plotlim,plotlim))
        ax[1].set(xlabel='x (Mpc)',ylabel='z (Mpc)',xlim=(-plotlim,plotlim),ylim=(-plotlim,plotlim))
        ax[2].set(xlabel='y (Mpc)',ylabel='z (Mpc)',xlim=(-plotlim,plotlim),ylim=(-plotlim,plotlim))

        ax[0].axis('square',fontsize=13)
        ax[1].axis('square')
        ax[2].axis('square')

        # setting labels
        labels = ['M31','MW','MW/M31 COM','M33 Orbit','M33','Stream Particle']

        #Plotting x,y plane
        ax[0].plot(0,0,'bo',markersize=3)
        ax[0].plot(MW_posx,MW_posy,'ro',markersize=3)
        ax[0].plot(MWM31_com_posx,MWM31_com_posy,'x',color='purple',markersize=5)
        ax[0].plot(M33_all_posx,M33_all_posy,'k',linewidth=1)
        ax[0].plot(M33_curr_posx,M33_curr_posy,'ko',markersize=3)
        ax[0].plot(streamx,streamy,'go',markersize=.5)

        #Plotting x,z plane
        ax[1].plot(0,0,'bo',markersize=3)
        ax[1].plot(MW_posx,MW_posz,'ro',markersize=3)
        ax[1].plot(MWM31_com_posx,MWM31_com_posz,'x',color='purple',markersize=5)
        ax[1].plot(M33_all_posx,M33_all_posz,'k',linewidth=1)
        ax[1].plot(M33_curr_posx,M33_curr_posz,'ko',markersize=3)
        ax[1].plot(streamx,streamz,'go',markersize=.5)

        #Plotting y,z plane
        ax[2].plot(0,0,'bo',markersize=3)
        ax[2].plot(MW_posy,MW_posz,'ro',markersize=3)
        ax[2].plot(MWM31_com_posy,MWM31_com_posz,'x',color='purple',markersize=5)
        ax[2].plot(M33_all_posy,M33_all_posz,'k',linewidth=1)
        ax[2].plot(M33_curr_posy,M33_curr_posz,'ko',markersize=3)
        ax[2].plot(streamy,streamz,'go',markersize=.5)

        # setting legend
        fig.legend(labels,bbox_to_anchor=[0.5,.925],loc='upper center',ncol=6,prop={'size':12},markerscale=3)

        # Saving movie frames
        plt.savefig('{}/local_group_video/{}.png'.format(master_path,ilbl))
        plt.tight_layout()
        plt.pause(0.001)

def velocity_dispersion(stream_data,r_min,r_max,leading=True,n_shells=20):
    """Function that calculates the radial and tangential velocity disperstion of
    a single leading or trailing stream. It is designed to caclulate the velocity dispersion
    at snap ID 465. This is beacuse snap 465 clearly shows leading and trailing stellar streams
    and is therefore a good frame to conduct this analysis on.
    
    Inputs:
    stream_data: numpy array of stream data at a single snapshot
    r_min: minimum radius to include stream particles (wrt M33 COM)
    r_max: maximum radius to include stream particles (wrt M33 COM)
    leading: If true, it calculates the velocity dispersion of the leading stream
             if false, it calculates the dispersion of the trailing stream
    n_shells: Number of radiall shells to calculate the dispersion within. 
              Just think of this as the number of bins
              
    Outputs:
    Plot of radial and tangential velocity dispersion as a function of radiall
    distance from M33.
    """

    #Centering particle positions and velocities wrt M33
    stream_xcen = stream_data['x'] - M33_data['x'][93]
    stream_ycen = stream_data['y'] - M33_data['y'][93]
    stream_zcen = stream_data['z'] - M33_data['z'][93]

    stream_vxcen = stream_data['vx'] - M33_data['vx'][93]
    stream_vycen = stream_data['vy'] - M33_data['vy'][93]
    stream_vzcen = stream_data['vz'] - M33_data['vz'][93]

    # These parameters are specific to snap 465. They are intended to exclud the wrap around parts of the stream
    if leading  == True:
        index1 = ((stream_ycen>=0) & (stream_ycen<=150))#conditional to be applied to data
        index2 = (stream_xcen<=25)
    else: 
        index1 = (stream_ycen<=0)
        index2 = (stream_xcen>=-25)

    #applying conditionals
    stream_xcen = stream_xcen[index1 & index2]
    stream_ycen = stream_ycen[index1 & index2]
    stream_zcen = stream_zcen[index1 & index2]

    stream_vxcen = stream_vxcen[index1 & index2]
    stream_vycen = stream_vycen[index1 & index2]
    stream_vzcen = stream_vzcen[index1 & index2]

    #creating a more convenient array for later
    stream_data = np.column_stack((stream_xcen,stream_ycen,stream_zcen,stream_vxcen,stream_vycen,stream_vzcen))
    
    # Initializing step size, this is so n_shells acts as a binning parameter
    bin_size = (r_max-r_min)/n_shells
    shell_radii = np.arange(r_min,r_max,bin_size)

    # Calculating radial distance for all stream particles from M33 COM
    r_vec = np.sqrt(stream_xcen**2 + stream_ycen**2 + stream_zcen**2)

    # Initializing array for storing stream data
    disperions = np.zeros((len(shell_radii),2))

    # Looping through shells
    for i in range(len(shell_radii)-1):        
        inner_shell = shell_radii[i] # radius of inner shell for ith bin
        outer_shell = shell_radii[i+1] # radius of inner shell for ith bin

        # First, I find stream particles at radial distances larger than the inner radius
        inner_shell_cond = (r_vec>=inner_shell) 
        inner_shell_data = stream_data[inner_shell_cond]
        inner_shell_r_vec = r_vec[inner_shell_cond]

        # Then I find particles at small radial distances than the outer radius
        outer_shell_cond = (inner_shell_r_vec<=outer_shell)
        outer_shell_data = inner_shell_data[outer_shell_cond]

        shell_data = outer_shell_data 

        # creating an array for r vector of all particles
        r_ = np.column_stack((shell_data[:,0],shell_data[:,1],shell_data[:,2]))
        r_mag = np.sqrt(r_[:,0]**2 + r_[:,1]**2 + r_[:,2]**2) #finding magnitude of r vector
        r_hat = np.array([r_[i,:]/r_mag[i] for i in range(len(r_mag))]) #finding unit r vector of all particles

        # Creating array for velocity vector for all partciles
        v = np.column_stack((shell_data[:,3],shell_data[:,4],shell_data[:,5]))
        v_mag = np.sqrt(v[:,0]**2 + v[:,1]**2 + v[:,2]**2) #finding magnitude of v vector
        v_r = np.array([np.dot(r_hat[i,:],v[i,:]) for i in range(len(v_mag))]) # dotting v vector with rhat to get velocity in r direction
        v_t = np.sqrt(v_mag**2 - v_r**2) #subtracting v_r from vmag to get velocity in tangential direction
        disperions[i,0] = np.std(v_r)#saving data to array
        disperions[i,1] = np.std(v_t)#saving data to array
    disperions_data = np.column_stack((shell_radii, disperions)) #adding shell radii to array
    path = "{}dispersion_data.txt".format(master_path) #saving data
    np.savetxt(path, disperions_data, fmt = "%11.4f"*3, comments='#',
               header="{:>11s}{:>11s}{:>11s}"\
                      .format('radii', 'sigma_r', 'sigma_t'))
