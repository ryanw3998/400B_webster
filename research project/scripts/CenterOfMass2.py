

"""
Ryan Webster
ASTR 400B
Modified CenterOfMass2 script
"""

# import modules
import numpy as np
# astropy provides unit system for astronomical calculations
import astropy.units as u
# astropy also offer a module for showinng results in a table
import astropy.table as tbl
# import previous HW functions
from ReadFile import Read

"""This is a modified version of the CenterOfMass2 script we built
in class. It includes the functions MassEnclosed and MassEnclosedTotal
from MassProfile.py. This allows me to read in files only one time. If I
used MassProfile.py, I would have to read them in twice."""
 
class CenterOfMass:
    """ Hold the COM position & velocity of a galaxy at a given snapshot """

    def __init__(self, filename, ptype):
        """ Initializing class with relevant data """
    
        # read in the file                                                                                             
        self.time, self.total, self.data = Read(filename)
        #print(self.time)                                                                                             
        
        self.filename = filename

        #create an array to store indexes of particles of desired Ptype                                                
        self.index = np.where(self.data['type'] == ptype)

        # This is a modification I made to get the positions of all particles in a galaxy.
        self.m0 = self.data['m']
        self.x0 = self.data['x']
        self.y0 = self.data['y']
        self.z0 = self.data['z']
         
        # store the mass, positions, velocities of only the particles of the given type 
        self.m = self.m0[self.index]
        self.x = self.x0[self.index]
        self.y = self.y0[self.index]
        self.z = self.z0[self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

    def COMdefine(self,a,b,c,m):
        """Function to compute the center of mass position or velocity generically
        Input: 
        array (a,b,c) of positions or velocities and the mass

        Returns: 
        3 floats  (the center of mass coordinates)         
        """
        # note: since all particles have the same                                                                      
        # mass, when we consider only one type,                                                                             
        # the below is equivalently np.sum(x)/len(x)     
        #print(a,b,c,m)                                                              

        # xcomponent Center of mass                                                                                    
        Acom = np.sum(a*m)/np.sum(m)
        # ycomponent Center of mass                                                                                    
        Bcom = np.sum(b*m)/np.sum(m)
        # zcomponent                                                                                                   
        Ccom = np.sum(c*m)/np.sum(m)
        return Acom, Bcom, Ccom
    
    # Modified For Homework 6
    def COM_P(self, delta, VolDec):
        """Iteratively determine the COM position
        Input:                                                                                                           
        delta (tolerance)  
        VolDec (value with which to decrease RMAX)
        
        Returns: 
        One vector, with rows indicating:                                                                                                                                                                            
        3D coordinates of the center of mass position (kpc)    
        """
    
        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x,self.y,self.z,self.m)
        # compute the magnitude of the COM position vector.                                                            
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)                                                                                 

        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position                                                                          
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2.0 + yNew**2.0 +zNew**2.0)

        # find the max 3D distance of all particles from the guessed COM
        # will re-start at a reduced radius specified by input VolDec                                                          
        RMAX = max(RNEW)/VolDec
        
        # pick an initial estimate for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume.                                
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
          # select all particles within the reduced radius (starting from original x,y,z, m)                         
            index2 = np.where(RNEW < RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius                                                                      
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position                                                                          
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)                                                                                 
            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by the specified decrement again                                                                
            RMAX = RMAX/VolDec                                                                                      

          # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM                                                                                     
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                         
            COMP = [XCOM,YCOM,ZCOM]

        # return the COM positon vector
        # set the correct units using astropy                                                                      
        # round all values 
        return np.round(COMP,2)*u.kpc

    def COM_P_test(self, delta, VolDec):
        """Iteratively determine the COM position
        Input:                                                                                                           
        delta (tolerance)  
        VolDec (value with which to decrease RMAX)
        
        Returns: 
        One vector, with rows indicating:                                                                                                                                                                            
        3D coordinates of the center of mass position (kpc)    
        """
    
        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x0,self.y0,self.z0,self.m0)
        # compute the magnitude of the COM position vector.                                                            
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)
        # print('init R', RCOM)                                                                                        


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position                                                                          
        xNew = self.x0 - XCOM
        yNew = self.y0 - YCOM
        zNew = self.z0 - ZCOM
        RNEW = np.sqrt(xNew**2.0 + yNew**2.0 +zNew**2.0)

        # find the max 3D distance of all particles from the guessed COM
        # will re-start at a reduced radius specified by input VolDec                                                          
        RMAX = max(RNEW)/VolDec
        
        # pick an initial estimate for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume.                                
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
          # select all particles within the reduced radius (starting from original x,y,z, m)                         
            index2 = np.where(RNEW < RMAX)
            x2 = self.x0[index2]
            y2 = self.y0[index2]
            z2 = self.z0[index2]
            m2 = self.m0[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius                                                                      
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position                                                                          
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # check this                                                                                               
            # print ("DIFF", diff)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by the specified decrement again                                                                
            RMAX = RMAX/VolDec
            # check this.                                                                                              
            #print ("maxR", maxR)                                                                                      

          # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM                                                                                     
            xNew = self.x0 - XCOM2
            yNew = self.y0 - YCOM2
            zNew = self.z0 - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                         
            COMP = [XCOM,YCOM,ZCOM]

        # return the COM positon vector
        # set the correct units using astropy                                                                      
        # round all values 
        return np.round(COMP,2)*u.kpc
    
    

    def COM_V(self, COMX,COMY,COMZ):
        """Return the COM velocity based on the COM position
        Input: 
        X, Y, Z positions of the COM
        
        Returns:
        3D Vector of COM Velocities
        """

        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position                              
        xV = self.x[:]*u.kpc - COMX
        yV = self.y[:]*u.kpc - COMY
        zV = self.z[:]*u.kpc - COMZ
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius                                                
        indexV = np.where(RV < RVMAX)

        # determine the velocity and mass of those particles within the mas radius                                     
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]

        # compute the center of mass velocity using those particles                                                    
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew, mnew)

        # create a vector to store the COM velocity                                                                    
        # set the correct units usint astropy                                                                          
        # round all values                                                                                             
        COMV = [VXCOM,VYCOM,VZCOM] 

        # return the COM vector  
        # set the correct units using astropy                                                                          
        # round all values   
        return np.round(COMV,2)*u.km/u.s

    def MassEnclosed(self, ptype, R, VolDec):
        """Function that determines the MassEnclosed of particles of a given type
        
        Input: 
        ptype: Particle type,  1=Halo, 2=Disk, 3=Bulge
        R: An Array of Radii within which to compute the mass enclosed. 

        Return: 
        An array with the Mass enclosed (units of Msun)
        """
        
        # Determine the COM position using Disk Particles
        # Disk Particles afford the best centroiding.
        # Set Delta = whatever you determined to be a good value in Homework 4.
        # Store the COM position of the galaxy
        # Set Delta = whatever you determined to be a good value in Homework 4.
        GalCOMP = self.COM_P(0.1,VolDec)
        # print(GalCOMP)
            
        # create an array to store indexes of particles of desired Ptype                                                
        index = np.where(self.data['type'] == ptype)

        # Store positions of particles of given ptype from the COMP. 
        xG = self.x0[index] - GalCOMP[0].value
        yG = self.y0[index] - GalCOMP[1].value
        zG = self.z0[index] - GalCOMP[2].value
            
        # Compute the mag. of the 3D radius
        rG = np.sqrt(xG**2 + yG**2 + zG**2)
            
        # store mass of particles of a given ptype
        mG = self.m0[index]
            
        # Array to store enclosed mass as a function of the input radius array
        Menc = np.zeros(np.size(R))
        # set up a while loop that continues until the end of the input radius array
        for i in range(np.size(R)):
            # Only want particles within the given radius
            indexR = np.where(rG < R[i])
            Menc[i] = np.sum(mG[indexR])*1e10         
        
        # return the array of enclosed mass with appropriate units
        return Menc*u.Msun
        
    
    def MassEnclosedTotal(self, R, VolDec):
        """Determine the total mass of each galaxy point

        Input:
        R: Radius from center of galaxy to desired end point
        Voldec: Needed for MassEnclosed function

        Returns: 
        Mencl: Mass in units Msun.
        """
            
        # Sum up all the mass of each component.
        Menc = self.MassEnclosed(1,R,VolDec) + self.MassEnclosed(2,R,VolDec) + self.MassEnclosed(3,R,VolDec)
    
        # Recall that M33 only has 2 components!  No bulge
        if 'M33' in self.filename:
            Menc = self.MassEnclosed(1,R,VolDec) + self.MassEnclosed(2,R,VolDec)  
          
        return Menc

def JacobiRadius(r_sep,m_sat,m_host):
    """Function to compute the Jacobi Radius given a separation between
    a galaxy and its satellite and their respective masses. Jacobi radius
    is calculated as Rj = R(Msat/2*Mhost(<R))^(1/3). Note: this equation 
    assumes an isothermal sphere for dark matter distribution. So the 
    returned value of radius might be under estimated.

    Intputs:
    r_sep: Host galaxy, satellite separation (kpc)
    m_host: Mass of host galaxy (Msun)
    m_sat: Mass of satellite galaxy (Msun)
    
    Returns:
    r_jacobi: Jacobi radius in kpc
    """

    #calculating jacobi radius
    r_jacobi = r_sep*(m_sat/(2*(m_host)))**(1/3)

    return r_jacobi

def separation(pos1,pos2):
    """Function to calculate the separation from one galaxy to another.
    
    Inputs:
    pos1: Numpy array of object 1's x,y,z position (kpc)
    pos2: Numpy array of object 2's x,y,z position (kpc)
    
    Returns:
    sep: Separation between the two bodies (kpc)
    """

    x = pos1[0] - pos2[0] #subtracting x values
    y = pos1[1] - pos2[1] #subtracting y values
    z = pos1[2] - pos2[2] #subtracting z values

    sep = np.sqrt(x**2 + y**2 + z**2) #Calulating the magnitude of separation vector

    return sep


