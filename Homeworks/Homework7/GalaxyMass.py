import numpy as np
from ReadFile import Read
from astropy.table import Table
from astropy.io import ascii

#from ParticleProperties import ParticleInfo

def CompnentMass(filename,particle_type):
    #Calculating total mass of specified type of particle

    #Inputs:
    #filename = path to file data file to analyze
    #particle_type = particle number correspodning to particle of interest (1 for Halo, 2 for Disk, 3 for Bulge)

    #Returns:
    #total mass of desired particle population

    time, num_part, data = Read(filename) #reading in data

    particle = data[data['type'] == particle_type] #Finding particle based on inputs particle_type
    tot_mass = np.sum(particle['m'])/100 #summing up all mass elements of a certain particle type. Converting to units of 10^12 M_sun
    return(np.around(tot_mass,3)) #returning mass


filename_MW = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/MW_000.txt' #defining file name
comp_mass_MW = [] #Creating empty list to fill with halo, disk and bulge mass
for i in range(1,4): #This loop is for computing the total mass of each particle type 
    comp_mass = CompnentMass(filename_MW,i) #calculating total mass for each particle type
    comp_mass_MW.append(comp_mass) #appending each componet to list

#Same process as above, except for M31
filename_M31 = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/M31_000.txt'
comp_mass_M31 = []
for i in range(1,4):
    comp_mass = CompnentMass(filename_M31,i)
    comp_mass_M31.append(comp_mass)

#Same process as above, except for M33
filename_M33 = '/Users/Ryan/Desktop/School/ASTR400B/400B_webster/M33_000.txt'
comp_mass_M33 = []
for i in range(1,4):
    comp_mass = CompnentMass(filename_M33,i)
    comp_mass_M33.append(comp_mass)

#Creating astropy table with appropriate columns
t = Table(names=('Galaxy Name', 'Halo Mass (10^12 Msun)', 'Disk Mass (10^12 Msun)', 'Buldge Mass (10^12 Msun)', 'Total (10^12 Msun)','fbar'), dtype=('S', 'f4', 'f4', 'f4', 'f4', 'f4'))

#Calculating fbar ratio for MW by simply adding up disk and bulge mass, and dividing by total system mass
mw_fbar = np.around((comp_mass_MW[1]+comp_mass_MW[2])/np.sum(comp_mass_MW),3)

#Same as above, but for M31
m31_fbar = np.around((comp_mass_M31[1]+comp_mass_M31[2])/np.sum(comp_mass_M31),3)

#Same as above, but for M33
m33_fbar = np.around((comp_mass_M33[1]+comp_mass_M33[2])/np.sum(comp_mass_M33),3)

#Calculating total group halo mass by simply adding up halo masses of MW, M31 and M33
group_halo_mass = np.around(comp_mass_MW[0] + comp_mass_M31[0] + comp_mass_M33[0],3)

#Same as above, but for total group disk mass
group_disk_mass = np.around(comp_mass_MW[1] + comp_mass_M31[1] + comp_mass_M33[1],3)

#Same as above, but for total group bulge mass
group_bulge_mass = np.around(comp_mass_MW[2] + comp_mass_M31[2] + comp_mass_M33[2],3)

#Adding up all group masses to find total mass
total_group_mass = group_bulge_mass+group_disk_mass+group_halo_mass

#Calculating f_bar for local group
f_bar = np.around((group_disk_mass+group_bulge_mass)/total_group_mass,3)

#Adding data to astropy table
t.add_row(('MW',comp_mass_MW[0],comp_mass_MW[1],comp_mass_MW[2],np.sum(comp_mass_MW),mw_fbar))
t.add_row(('M31',comp_mass_M31[0],comp_mass_M31[1],comp_mass_M31[2],np.sum(comp_mass_M31),m31_fbar))
t.add_row(('M33',comp_mass_M33[0],comp_mass_M33[1],comp_mass_M33[2],np.sum(comp_mass_M33),m33_fbar))
t.add_row(('Group',group_halo_mass,group_disk_mass,group_bulge_mass,total_group_mass,f_bar))

#Saving table
ascii.write(t,'comp_masses', format='latex')