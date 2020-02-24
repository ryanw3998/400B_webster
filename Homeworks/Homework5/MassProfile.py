
#HOMEWORK #5 
#Ryan Webster
#Collaborators: Mackenzie James


#Importing modules
import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterOfMass import CenterOfMass
import matplotlib.pyplot as plt
from astropy.constants import G 
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

class MassProfile:
#Class for determining the mass profile of a galaxy

    def __init__(self, galaxy, Snap):
        #Initializing class with a galaxy name and snap number

        ilbl = '000' + str(Snap)
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        #The code above lets us simply input a galaxies name, rather than
        #a file name when using the class

        self.gname = galaxy # storing the galaxy's name

        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)                                                                                             
    
    def MassEnclosed(self, ptype, radii):
        #Function for finding the mass enclosed within an array of radii

        #Inputs:
            #ptype: The particle type 
            #radii: A numpy array of radii
        #Returns:
            #An array of mass enclosed within each given radii, in units of 
            # 1e10 Msun

        COM = CenterOfMass(self.filename,2) #Finding the center of mass
        COMP = COM.COM_P(0.1)

        self.index = np.where(self.data['type'] == ptype) #Selecting particle type
        self.m = self.data['m'][self.index] #getting particle masses
        self.x = self.data['x'][self.index]*u.kpc #getting position information
        self.y = self.data['y'][self.index]*u.kpc
        self.z = self.data['z'][self.index]*u.kpc

        x_coord = self.x - COMP[0] #Calculating positional coordinates of each particle
        y_coord = self.y - COMP[1]
        z_coord = self.z - COMP[2]

        r_coord = np.sqrt(x_coord**2+y_coord**2+z_coord**2) #Finding magnitude of position

        masses = np.zeros(len(radii)) #Initializing array

        for i in range(len(radii)):#For each radius in radii,
            encl = np.where(r_coord <= radii[i])  #finding particles within that radius
            enclosedmass = self.m[encl] #finding their mass
            masses[i] = np.sum(enclosedmass) #summing up their masses

        return masses*1e10*u.Msun

    def MassEnclosedTotal(self,radii):
        #Function that returns the total mass of all particle types within
        #each radius in radii

        #Inputs:
            #radii: A numpy array of radii

        #Returns:
            #Mass of all particle types within each radius in radii

        if self.gname == 'M33': #This conditional is because M33 doesn't have a bulge
            halo = self.MassEnclosed(1,radii)#Calling upon previously written MassEnclosed function
            disk = self.MassEnclosed(2,radii)
            return(halo+disk) #Summing up mass componets
        else:
            halo = self.MassEnclosed(1,radii)
            disk = self.MassEnclosed(2,radii)
            bulge = self.MassEnclosed(3,radii)
            return(halo+disk+bulge) #Summing up mass componets

    def HernquistMass(self, radius, a, Mhalo):
        #Function for calculating the mass enclosed within a certain radius using the 
        #theoretical Hernquist Mass profile

        #Inputs:
            #radius: Radius for enclosed mass to be calculated
            #a: Scale factor
            #Mhalo: Halo mass
        #Returns:
            #Mass enclosed within given radius based on Hernquist Mass Profile
            #in units of 1e12 Msun

        M = np.round(Mhalo*radius**2 / (a + radius)**2,2) #calculating mass
        #This function uses M(r) = (M_halo*r^2)/(a+r)^2 to calculate mass profile

        return M*u.Msun

    def CircularVelocity(self, ptype, radii):
        #Function for calculating circular velocity of particles with given radii using
        #v = sqrt(GM/R)

        #Inputs:
            #ptype: Particle type to calculate the circular velocity of
            #radii: Numpy array of radii to calculate the circular velocity of 

        #Returns:
            #circular velocity of particle in km/s

        masses = self.MassEnclosed(ptype,radii) #using mass enclosed to get mass for velocity eq
        circvel = np.round(np.sqrt(G*masses/radii),2) #calculating circular velocity

        return circvel*u.km/u.s

    def CircularVelocityTotal(self, radii):
        #Function for calculating the circular velocity of all particle types of a galaxy
        #Uses v =sqrt(GM/R)

        #Inputs:
            #radii: Numpy array of radii
        #Returns:
            #circular velocity of all particle types in km/s
        if self.gname == 'M33': #Conditional for M33, which doesn't have a buldge
            halo = self.MassEnclosed(1,radii)
            disk = self.MassEnclosed(2,radii)
            totalmass = halo+disk #Adding together each componet mass
        else:
            halo = self.MassEnclosed(1,radii) #Extracting each componet mass
            disk = self.MassEnclosed(2,radii)
            buldge = self.MassEnclosed(3, radii)
            totalmass = halo+disk+buldge #Adding together each componet mass

        circveltotal = np.round(np.sqrt(G*totalmass/radii),2) #calculating circular velocity

        return circveltotal*u.km/u.s

    def HernquistVCirc(self, radii, a, Mhalo):
        #Function for calculating the cirular speed based off the 
        #Hernquist mass profile. Uses v = sqrt(GM/R)

        #Inputs:
            #radii: Numpy array of radii
            #a: scale factor for Hernquist Mass Profile
            #Mhalo: For Hernquist mass profile
        #Outputs:
            #circular velocities for each radii

        hernquist_mass = self.HernquistMass(radii, a, Mhalo)#Calling Hernquist mass profile from earlier
        circvel = np.sqrt(G*hernquist_mass/radii) #calculating circular velocity

        return np.round(circvel,2) *u.km/u.s

#******** Section 8 ********

# Data for Milky Way
r = np.arange(0.25, 30.5, 1.5)*u.kpc
MWmp = MassProfile('MW',0)
MW_mass_halo = MWmp.MassEnclosed(1,r)
MW_mass_disk = MWmp.MassEnclosed(2,r)
MW_mass_bulge = MWmp.MassEnclosed(3,r)
MW_mass_total = MWmp.MassEnclosedTotal(r)
MW_mhalo = 1.975*1e12 *u.Msun
MW_a = 65 *u.kpc
MW_hernmass = MWmp.HernquistMass(r, MW_a, MW_mhalo)

# Plotting mass profile for Milky Way
plt.figure(1)
plt.title('Section 8: MW Mass Profile')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Mass Enclosed ($M_{sun}$)')
plt.semilogy()
plt.plot(r,MW_mass_halo,'--',label='Halo')
plt.plot(r,MW_mass_disk,'--',label='Disk')
plt.plot(r,MW_mass_bulge,'--',label='Bulge')
plt.plot(r,MW_mass_total,'--',label='Total')
plt.legend()
#plt.savefig('Plots/MW_mass_profile.png',dpi=300)
plt.show()

# Plotting Milky Way vs Hernquist Mass Profile

plt.figure(2)
plt.title('Section 8: MW Halo vs Hernquist')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Mass Enclosed ($M_{sun}$)')
plt.semilogy()
plt.plot(r,MW_mass_halo,'--',label='Halo')
plt.plot(r,MW_hernmass,'--',label='Hernquist Profile, a={}'.format(MW_a))
plt.legend()
#plt.savefig('Plots/MW_hernquist.png',dpi=300)
plt.show()

# Data for M31 mass profile

r = np.arange(0.25, 30.5, 1.5)*u.kpc
M31mp = MassProfile('M31',0)
M31_mass_halo = M31mp.MassEnclosed(1,r)
M31_mass_disk = M31mp.MassEnclosed(2,r)
M31_mass_bulge = M31mp.MassEnclosed(3,r)
M31_mass_total = M31mp.MassEnclosedTotal(r)
M31_mhalo = 1.921*1e12 * u.Msun
M31_a = 60 *u.kpc
M31_hernmass = M31mp.HernquistMass(r, M31_a, M31_mhalo)

# Plotting M31 mass profile 

plt.figure(3)
plt.title('Section 8: M31 Mass Profile')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Mass Enclosed ($M_{sun}$)')
plt.semilogy()
plt.plot(r,M31_mass_halo,'--',label='Halo')
plt.plot(r,M31_mass_disk,'--',label='Disk')
plt.plot(r,M31_mass_bulge,'--',label='Bulge')
plt.plot(r,M31_mass_total,'--',label='Total')
plt.legend()
#plt.savefig('Plots/M31_mass_profile.png',dpi=300)
plt.show()

# Plotting M31 vs Hernquist mass profile

plt.figure(4)
plt.title('Section 8: M31 Halo vs Hernquist')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Mass Enclosed ($M_{sun}$)')
plt.semilogy()
plt.plot(r,M31_mass_halo,'--',label='Halo')
plt.plot(r,M31_hernmass,'--',label='Hernquist Profile, a={}'.format(M31_a))
plt.legend()
#plt.savefig('Plots/M31_hernquist.png',dpi=300)
plt.show()

# Data for M33 mass profile

r = np.arange(0.25, 30.5, 1.5)*u.kpc
M33mp = MassProfile('M33',0)
M33_mass_halo = M33mp.MassEnclosed(1,r)
M33_mass_disk = M33mp.MassEnclosed(2,r)
M33_mass_bulge = M33mp.MassEnclosed(3,r)
M33_mass_total = M33mp.MassEnclosedTotal(r)
M33_mhalo = 1.921*1e12 * u.Msun
M33_a = 90 *u.kpc
M33_hernmass = M33mp.HernquistMass(r, M33_a, M33_mhalo)

#Plotting M33 mass profile

plt.figure(5)
plt.title('Section 8: M33 Mass Profile')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Mass Enclosed ($M_{sun}$)')
plt.semilogy()
plt.plot(r,M33_mass_halo,'--',label='Halo')
plt.plot(r,M33_mass_disk,'--',label='Disk')
plt.plot(r,M33_mass_bulge,'--',label='Bulge')
plt.plot(r,M33_mass_total,'--',label='Total')
plt.legend()
#plt.savefig('Plots/M33_mass_profile.png',dpi=300)
plt.show()

#Plotting M33 vs Hernquist mass profile

plt.figure(6)
plt.title('Section 8: M33 Halo vs Hernquist')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Mass Enclosed ($M_{sun}$)')
plt.semilogy()
plt.plot(r,M33_mass_halo,'--',label='Halo')
plt.plot(r,M33_hernmass,'--',label='Hernquist Profile, a={}'.format(M33_a))
plt.legend()
#plt.savefig('Plots/M33_hernquist.png',dpi=300)
plt.show()

#******** Section 9 ********

#Data for Milky Way rotation curve
r = np.arange(0.25, 30.5, 1.5)*u.kpc
MWmp = MassProfile('MW',0)
MW_vcir_halo = MWmp.CircularVelocity(1,r)
MW_vcir_disk = MWmp.CircularVelocity(2,r)
MW_vcir_bulge = MWmp.CircularVelocity(3,r)
MW_vcir_total = MWmp.CircularVelocityTotal(r)
MW_mhalo = 1.975*1e12 *u.Msun
MW_a = 65 *u.kpc
MW_herncirv = MWmp.HernquistVCirc(r, MW_a, MW_mhalo)

#Plotting MW rotation curve 

plt.figure(7)
plt.title('Section 9: MW Rotaion Curve')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Circular Velocity ($km/s$)')
plt.semilogy()
plt.plot(r,MW_vcir_halo,'--',label='Halo')
plt.plot(r,MW_vcir_disk,'--',label='Disk')
plt.plot(r,MW_vcir_bulge,'--',label='Bulge')
plt.plot(r,MW_vcir_total,'--',label='Total')
plt.legend()
#plt.savefig('Plots/MW_rotation_curve.png',dpi=300)
plt.show()

#Plotting Milky Way halo vs Hernquist rotation curve

plt.figure(8)
plt.title('Section 9: MW Halo vs Hernquist')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Circular Velocity ($km/s$)')
plt.semilogy()
plt.plot(r,MW_vcir_halo,'--',label='Halo')
plt.plot(r,MW_herncirv,'--',label='Hernquist Profile, a={}'.format(MW_a))
plt.legend()
#plt.savefig('Plots/MW_hernquist_rot_curve.png',dpi=300)
plt.show()

#Data for M31 rotation curve 

r = np.arange(0.25, 30.5, 1.5)*u.kpc
M31mp = MassProfile('M31',0)
M31_vcir_halo = M31mp.CircularVelocity(1,r)
M31_vcir_disk = M31mp.CircularVelocity(2,r)
M31_vcir_bulge = M31mp.CircularVelocity(3,r)
M31_vcir_total = M31mp.CircularVelocityTotal(r)
M31_mhalo = 1.975*1e12 *u.Msun
M31_a = 60 *u.kpc
M31_herncirv = M31mp.HernquistVCirc(r, M31_a, M31_mhalo)

#Plotting M31 rotation curve

plt.figure(9)
plt.title('Section 9: M31 Rotation Curve')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Circular Velocity ($km/s$)')
plt.semilogy()
plt.plot(r,M31_vcir_halo,'--',label='Halo')
plt.plot(r,M31_vcir_disk,'--',label='Disk')
plt.plot(r,M31_vcir_bulge,'--',label='Bulge')
plt.plot(r,M31_vcir_total,'--',label='Total')
plt.legend()
#plt.savefig('Plots/M31_rotation_curve.png',dpi=300)
plt.show()

#Plotting M31 halo vs Hernquist rotation curve

plt.figure(10)
plt.title('Section 9: M31 Halo vs Hernquist')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Circular Velocity ($km/s$)')
plt.semilogy()
plt.plot(r,M31_vcir_halo,'--',label='Total')
plt.plot(r,M31_herncirv,'--',label='Hernquist Profile, a={}'.format(M31_a))
plt.legend()
#plt.savefig('Plots/M31_hernquist_rot_curve.png',dpi=300)
plt.show()

#Data for M33 rotation curve

r = np.arange(0.25, 30.5, 1.5)*u.kpc
M33mp = MassProfile('M33',0)
M33_vcir_halo = M33mp.CircularVelocity(1,r)
M33_vcir_disk = M33mp.CircularVelocity(2,r)
M33_vcir_total = M33mp.CircularVelocityTotal(r)
M33_mhalo = 1.975*1e12 *u.Msun
M33_a = 90 *u.kpc
M33_herncirv = M33mp.HernquistVCirc(r, M33_a, M33_mhalo)

#Plotting M33 rotation curve

plt.figure(11)
plt.title('Section 9: M33 Rotation Curve')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Circular Velocity ($km/s$)')
plt.semilogy()
plt.plot(r,M33_vcir_halo,'--',label='Halo')
plt.plot(r,M33_vcir_disk,'--',label='Disk')
plt.plot(r,M33_vcir_total,'--',label='Total')
plt.legend()
#plt.savefig('Plots/M33_rotation_curve.png',dpi=300)
plt.show()

#Plotting M33 halo vs Hernquist rotation curve

plt.figure(12)
plt.title('Section 9: M33 Halo vs Hernquist')
plt.xlabel('Radius ($kpc$)')
plt.ylabel('Circular Velocity ($km/s$)')
plt.semilogy()
plt.plot(r,M33_vcir_halo,'--',label='Total')
plt.plot(r,M33_herncirv,'--',label='Hernquist Profile, a={}'.format(M33_a))
plt.legend()
#plt.savefig('Plots/M33_hernquist_rot_curve.png',dpi=300)
plt.show()