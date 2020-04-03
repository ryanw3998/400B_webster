
# Homework 7 code
# Ryan Webster
# Collaborators: Jimmy Lilly, Madison Walder, Mackenzie James, Sean Cunningham, and Apollo (hes my cat)


import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from IPython.display import Latex
from CenterOfMass2 import CenterOfMass

# # M33AnalyticOrbit
class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        """Initializing class with a filename
        The filename is the path to where you want to store the analytical data"""
        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
       
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33COM = CenterOfMass("M33_000.txt", 2)
        # **** store the position VECTOR of the M33 COM 
        posM33 = M33COM.COM_P(0.1,4) #I use .value later, was getting error in COM_V code
        # **** store the velocity VECTOR of the M33 COM 
        velM33 = M33COM.COM_V(posM33[0],posM33[1],posM33[2]) 
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31COM = CenterOfMass('M31_000.txt', 2)
        # **** store the position VECTOR of the M31 COM 
        posM31 = M31COM.COM_P(0.1,2) #I use .value later, was getting error in COM_V code
        # **** store the velocity VECTOR of the M31 COM 
        velM31 = M31COM.COM_V(posM31[0],posM31[1],posM31[2])
        
        # relative position and velocity VECTORS of M33
        self.r0 = posM33.value - posM31.value
        self.v0 = velM33.value - velM31.value
        
        ### get the mass of each component in M31 
        ### disk
        self.rdisk = 5. #kpc
        self.Mdisk = 0.12*10**12 #Msun
        ### bulge
        self.rbulge = 1. #kpc
        self.Mbulge = 0.019*10**12 #Msun
        # Halo
        self.rhalo = 60. #kpc
        self.Mhalo = 1.921*10**12 #Msun
    
    def HernquistAccel(self, M, ra, r): 
        """ Function to calculate the hernquist acceleration of the halo and bulge
        Inputs:
        M = mass of halo or bulge in Msun (but without Astropy units)
        ra = scale legth of halo or bulge in kpc (but without Astropy units)
        r = position vector
        Outputs:
        Acceleration vector, as calculated by Hernquist equation """
        
        ###  Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        
        ###  Store the Acceleration
        Hern =  -self.G*M/(rmag *(ra + rmag)**2) * r 
        
        return Hern
    
    
    def MiyamotoNagaiAccel(self, M, rd, r):
        """ Function to calculate the Miyamoto-Nagai acceleration of the disk
        Inputs:
        M = mass of disk in Msun (but without Astropy units)
        ra = scale legth of halo or bulge in kpc (but without Astropy units)
        r = position vector 
        Outputs:
        Acceleration vector, as calculated by Miyamoto-Nagai equation """
        R = np.sqrt(r[0]**2 + r[1]**2) #Finding magnitude of x and y compnets
        zd = rd/5. #Calculating "zd"
        B = rd + np.sqrt(r[2]**2 + zd**2) #Calclating "B"
        zstuff = 1/np.sqrt(r[2]**2 + zd**2) #Calculating stuff that only appears in z componet
        MNa = -self.G*M/(R**2+B**2)**1.5 * r * np.array([1,1,zstuff]) #Putting it all together

        return MNa
     
    
    def M31Accel(self,r): 
        """ Function to sum componets of halo, bulge, and disk accelerations
        Inputs:
        r = position vector
        
        Outputs:
        Summed acceleration vector """

        ### Call the previous functions for the halo, bulge and disk
        halo_acc = self.HernquistAccel(self.Mhalo,self.rhalo,r) #Calculating halo acceleration
        bulge_acc = self.HernquistAccel(self.Mbulge,self.rbulge,r) #""" bulge """
        disk_acc = self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,r)#""" disk """
        sum1 = np.add(halo_acc,bulge_acc) #summing halo and bulge first
        acc_sum = np.add(sum1,disk_acc) #then adding disk acceleration

        return acc_sum
    
    
    
    def LeapFrog(self,r,v,dt):
        """ Funtion that impliments leap frog numerical integration method
        Inputs:
        r = position vector
        v = velocity vector
        dt = timestep
        
        Ouputs:
        Next interval in integration method for position and velocity vectors """

        rhalf = r + np.asarray(v)*(dt/2) #Taking a half step forward with positional vector
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + self.M31Accel(rhalf)*dt
        # predict the final position using the average of the current velocity and the final velocity
        rnew = r + 0.5*(v+vnew)*dt
        
        return rnew,vnew
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ Integrating the orbit using leap frog algorithm
        Inputs:
        t0 = starting time
        dt = timestep to use for integration
        tmax = ending time
        
        Outputs:
        Text file contaning the position and velocity vector at each timestep 
        from t0 to tmax"""

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size
        orbit = np.zeros((int(tmax/dt)+2,7))
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        
        # initialize a counter for the orbit.  
        i = 1 
        
        # start the integration 
        while (t<tmax): 
            # **** advance the time by one timestep, dt
            t += dt
            orbit[i-1,0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            r = [orbit[i-1,1], orbit[i-1,2], orbit[i-1,3]]
            v = [orbit[i-1,4], orbit[i-1,5], orbit[i-1,6]]
            rnew, vnew = self.LeapFrog(r,v,dt)
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i] = t, * tuple(rnew), *tuple(vnew)
            # **** update counter i
            i += 1
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


#****************************************
#************* Questions ****************

#Question 1: Comparing Analytical Separation to Simulation
#Creating instance of OrbitIntegration
file = 'OrbitIntergartion.txt'
orbit = M33AnalyticOrbit(file)
orbit.OrbitIntegration(0,0.1,10)

#Reading in analytical data
data = np.loadtxt(file)
last = len(data)-1 #Getting rid of last row in data, it was all zeros
#Calculating magnitude of distance vector between M33 and M31
mag_r = np.sqrt(data[:last,1]**2+data[:last,2]**2,data[:last,3]**2)
mag_v = np.sqrt(data[:last,4]**2+data[:last,5]**2,data[:last,6]**2)

#Reading in Orbit data from assignment 6
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

#Using mag_diff for M33 and M31
M33_M31_r,M33_M31_vr = mag_diff(M33_data,M31_data)
#Plotting simulation and analytical orbit against each other
plt.figure(1)
plt.title("M31 M33 Separation (0-10 Gyr)")
plt.xlabel('Time (Gyr)')
plt.ylabel('Separation (kpc)')
plt.plot(data[:last,0],mag_r,label='Analytical')
plt.plot(M31_data['t'],M33_M31_r,label='Simulated')
plt.legend()
plt.savefig('M33_M31_separation.png',dpi=350)

plt.show()

plt.figure(2)
plt.title("M31 M33 Velocity (0-10 Gyr)")
plt.xlabel('Time (Gyr)')
plt.ylabel('Velocity (km/s)')
plt.plot(data[:last,0],mag_v,label='Analytical')
plt.plot(M31_data['t'],M33_M31_vr,label='Simulated')
plt.legend()
plt.savefig('M33_M31_velocity.png',dpi=350)

plt.show()

#Question 2: 
"""The plots look very different from each other.
The separation simulation shows the two galaxies 
ocsillating between ~100 and ~200 kpc and slowly getting
closer to each other as the simulation goes on. The analytical 
plot shows the two starting at ~150 kpc, getting to about ~100 kpc
and then the separation increasing dramatically to ~400 kpc before
starting to decrease again. A similar story is seen is the velocity
plot. To me, the trajectory of the analytical plot indicates that 
M33 got very close to M31, and was then flung out to a large radius, 
similar to a comet on a highly eccentirc orbit around the sun."""

#Question 3:
"""The anaylitical approach treated M31's disk, bulge, and halo 
componets individually, whereas the simulation broke each galaxy 
componet down in much smaller chuncks. This means the analytical
approach is much more general in its approach, and not as 
accurate. This could account for some of the discrepancies between
methods. """

#Question 4:
"""The other factor at play is the lack of MW in the system. Being
the second most massive body in the local group, MW will play a
significant role in the dynamics of the M31-M33 interaction. Not 
having the MW is likly the biggest reason for the discrepancy between
the simulation and analytical approach. To include its effects,
this approach above would need to include the interatction between 
MW and M33. 
"""