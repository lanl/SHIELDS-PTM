"""
This file contains routines necessary to create input files for the ptm simulation. See in-line comments
for the description of various parameters.

Usage:

    The simplest case is used to generate a single particle drift trace for a pre-defined test:

    import ptm_input
    p=ptm_input.ptm_input_creator()
    p.create_input_files()

This is an object-oriented version of a procedural code of the same name.

Jesse Woodroffe
5/30/2017

"""

import numpy as np
import os

class ptm_input_creator(object):
    """
    An object for creating input data files for the SHIELDS-PTM code framework
    """

    __ckm = 2.99792458e5         # Speed of light in km/s

    # These are internal parameter sets which differ based on certain arguments.

    __pkeys=['runid','nparticles','ndim','nx','ny','nz','itrace','ifirst','ilast','ntot','dtin','dtout',
             'istep','iswitch','iphase','nphase','charge','mass','tlo','thi','itraj']

    __dkeys=[['idens','x0','y0','z0'],
             ['idens','xmin','xmax','ymin','ymax','zmin','zmax'],
             ['idens','r0','mltmin','mltmax']]

    __vkeys=[['idist','ekev','alpha','phi'],
             ['idist','vtperp','vtpara','phi'],
             ['idist','nenergy','npitch','phi','emin','emax','pamin','pamax','xsource'],
             ['idist','nenergy','npitch','phi','xsource']]

    def __init__(self,runid=1,idensity=1,ivelocity=1):
        """
        Initialize to default values for specified parameter sets
        """

        self.__pdict={}
        self.__ddict={}
        self.__vdict={}

        self.__pdict['runid']=runid
        self.__pdict['nparticles']=1   # Number of particles to use in simulation
        self.__pdict['ndim']=3         # Number of spatial dimensions of input data (2 or 3)
        self.__pdict['nx']=75          # Number of points in X grid
        self.__pdict['ny']=75          # Number of points in Y grid
        self.__pdict['nz']=75          # Number of points in Z grid
        self.__pdict['itrace']=-1      # 1 = forward in time, -1 = backwards in time
        self.__pdict['ifirst']=1       # Timestep of first fields snapshot: e.g., start at bx3d_0005.bin if ifirst=5
        self.__pdict['ilast']=2        # Timestep of last fields snapshot: e.g., finish at bx3d_0009.bin if ilast=9
        self.__pdict['ntot']=2         # Total number of time steps in tgrid file
        self.__pdict['dtin']=3600.0    # Time between snapshots (seconds). If dtin<0, read times from tgrid.bin
        self.__pdict['dtout']=1.0      # Time between data outputs (seconds)
        self.__pdict['istep']=1        # Integrator to use: 1=RK4, 2=Adaptive RKSuite
        self.__pdict['iswitch']=-1     # Particle equations to solve: -1=guiding center, 0=dynamic switching, 1=full orbit
        self.__pdict['iphase']=3       # Gyrophase switching method: 1=random, 2=brute force search, 3=gradient method
        self.__pdict['nphase']=32      # Number of points for brute force search (minimize difference in B)
        self.__pdict['charge']=-1.0    # Particle charge in multiples of fundamental
        self.__pdict['mass']=1.0       # Particle mass in multiples of the electron mass
        self.__pdict['tlo']=0.0        # Lowest time of particle trace
        self.__pdict['thi']=3600.0     # Highest time of particle trace
        self.__pdict['itraj']=0        # Write particle trajectories if in flux map mode: 0 = no, 1 = yes

        self.__vdict['idist']=ivelocity

        if(ivelocity==1):
            # Ring
            self.__vdict['idist']=1
            self.__vdict['ekev']=100.0     # Particle energy in keV
            self.__vdict['alpha']=90.0     # Particle pitch angle in degrees
            self.__vdict['phi']=180.0      # Particle gyrophase angle in degrees; set negative to seed randomly
        elif(ivelocity==2):
            # Bi-Maxwellian
            self.__vdict['idist']=2
            self.__vdict['vtperp']=0.25*__ckm  # Perpendicular thermal velocity in km/s
            self.__vdict['vtpara']=0.25*__ckm  # Parallel thermal velocity in km/s
            self.__vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
        elif(ivelocity==3):
        # Uniform flux map mode
            self.__vdict['idist']=3
            self.__vdict['nenergy']=34       # Number of energies to use in discrete flux map
            self.__vdict['npitch']=31        # Number of pitch angles to use in discrete flux map
            self.__vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
            self.__vdict['emin']=1.0         # Lowest energy of flux map (keV)
            self.__vdict['emax']=100.0       # Highest energy of flux map (keV)
            self.__vdict['pamin']=0.0        # Lowest pitch angle of flux map (degrees)
            self.__vdict['pamax']=90.0       # Highest pitch angle of flux map (degrees)
            self.__vdict['xsource']=-12.0    # XGSM distance where integration should be terminated (RE)
        elif(ivelocity==4):
        # User-specified flux map mode
            self.__vdict['idist']=3
            self.__vdict['nenergy']=34       # Number of energies to use in discrete flux map
            self.__vdict['npitch']=31        # Number of pitch angles to use in discrete flux map
            self.__vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
            self.__vdict['xsource']=-12.0    # XGSM distance where integration should be terminated (RE)
        else:
            raise Exception('Option ivelocity = '+str(ivelocity)+' is not supported.')

        if(idensity==1):
            # Single point
            self.__ddict['idens']=1
            self.__ddict['x0']=-5.0      # Initial X-position (RE)
            self.__ddict['y0']=0.0       # Initial Y-position (RE)
            self.__ddict['z0']=0.0       # Initial Z-position (RE)
        elif(idensity==2):
            # Cubic region
            self.__ddict['idens']=2
            self.__ddict['xmin']=-7.0    # Minimum X-boundary of seeding region (RE)
            self.__ddict['xmax']=-6.0    # Maximum X-boundary of seeding region (RE)
            self.__ddict['ymin']=-0.5    # Minimum Y-boundary of seeding region (RE)
            self.__ddict['ymax']= 0.5    # Maximum Y-boundary of seeding region (RE)
            self.__ddict['zmin']=-0.5    # Minimum Z-boundary of seeding region (RE)
            self.__ddict['zmax']=0.5     # Maximum Z-boundary of seeding region (RE)
        elif(idensity==3):
            # Radial distance
            self.__ddict['idens']=3
            self.__ddict['r0']=6.6       # Radial distance for seeding (RE)
            self.__ddict['mltmin']=-6.0  # Minimum local time for seeding region (HOURS)
            self.__ddict['mltmax']=6.0   # Maximum local time for seeding region (HOURS)

    def get_interactive_input(self):
        """
        Set parameter inputs using interactive prompts. These quantities need to have more informative messages
        associated with them, but the infrastructure is functional.
        """

        prompt1="Please follow prompts to set parameters for PTM simulation"
        prompt2="Step 1: Set global simulation parameters"
        prompt3="Step 2: Set parameters for configuration space"
        prompt4="Step 3: Set parameters for velocity space"
        prompt5="Create input files? (Y/N)\n"

        print('*'*len(prompt1))
        print(prompt1)
        print('*'*len(prompt1))

        runid=input("\nEnter run identification number (positive integer)\n")
        idens=input("\nEnter configuration type (1=point,2=cubic domain,3=cirular region)\n")
        idist=input("\nEnter distribution type (1=ring,2=Maxwellian,3=Uniform FluxMap,4=Non-uniform FluxMap)\n")

        # Get parameter lists for these options
        self.__init__(runid=runid,idensity=idens,ivelocity=idist)

        # Set parameters for global simulation
        print('*'*len(prompt2))
        print(prompt2)
        print('*'*len(prompt2))
        for key in np.setxor1d(['runid'],self.__pkeys):
            self.__pdict[key]=input("\n"+key+"=\n")

        # Set configuration space parameters
        print('\n'+'*'*len(prompt3))
        print(prompt3)
        print('*'*len(prompt3))
        for key in np.setxor1d(['idens'],self.__dkeys[self.__ddict['idens']-1]):
            self.__ddict[key]=input("\n"+key+"=\n")

        # Set velocity space parameters
        print('\n'+'*'*len(prompt4))
        print(prompt4)
        print('*'*len(prompt4))
        for key in np.setxor1d(['idist'],self.__vkeys[self.__vdict['idist']-1]):
            self.__vdict[key]=input("\n"+key+"=\n")

        # Decide whether to create output files

        print('\n'+'*'*len(prompt5))
        yesno = input(prompt5)

        if(yesno[0].upper()=='Y'):
            self.create_input_files()

        return

    def set_parameters(self,**kwargs):
        """
        Change internal parameter settings
        """

        if('runid' in kwargs.keys()):
            self.__pdict['runid']=kwargs['runid']

        pdict=self.__pdict.copy()
        ddict=self.__ddict.copy()
        vdict=self.__vdict.copy()

        replace_D = False
        replace_V = False

        if('idens' in kwargs.keys()):
            idens=kwargs['idens']
            replace_D = True
        else:
            idens=ddict['idens']

        if('idist' in kwargs.keys()):
            idist=kwargs['idist']
            replace_V = True
        else:
            idist=vdict['idist']

        # Change parameter sets if appropriate and put back values that should not change
        if(replace_D or replace_V):

            self.__init__(runid=pdict['runid'],idensity=idens,ivelocity=idist)

            for key,value in pdict.items():
                if(key in self.__pkeys):
                    self.__pdict[key]=value
            for key,value in ddict.items():
                if(key in self.__dkeys):
                    self.__ddict[key]=value
            for key,value in vdict.items():
                if(key in self.__vkeys):
                    self.__vdict[key]=value

        # Set values from user input
        for key,value in kwargs.items():
            if(key in self.__pkeys):
                self.__pdict[key]=value
            for dkeys in self.__dkeys:
                if(key in dkeys):
                    self.__ddict[key]=value
            for vkeys in self.__vkeys:
                if(key in vkeys):
                    self.__vdict[key]=value

        return

    def print_settings(self):
        """
        Print out a summary of internal parameters
        """

        print("*********************")
        print("Simulation Parameters")
        print("*********************")
        for key in self.__pkeys:
            print('{:<12}{:>8}'.format(key, self.__pdict[key]))
        print("")
        print("******************************")
        print("Configuration Space Parameters")
        print("******************************")
        for key in self.__dkeys[self.__ddict['idens']-1]:
            print('{:<12}{:>8}'.format(key, self.__ddict[key]))
        print("")
        print("*************************")
        print("Velocity Space Parameters")
        print("*************************")
        for key in self.__vkeys[self.__vdict['idist']-1]:
            print('{:<12}{:>8}'.format(key, self.__vdict[key]))

        return

    def create_input_files(self,filedir=''):
        """
        Write the output files
        """

        # Handle a couple of special cases
        if(len(filedir)>0):
            if(filedir=='.'):
                # User specified present directory using dot; strip the dot
                fileir=''
            elif(filedir[-1]!='/'):
                # User forgot to add a slash at the end of the directory name; add a slash
                filedir+='/'

        # Output filenames
        #TODO:check whether ptm_input subdirectory needs to be specified here
        pfname='{}ptm_input/ptm_parameters_{:04d}.txt'.format(filedir,self.__pdict['runid'])
        dfname='{}ptm_input/dist_density_{:04d}.txt'.format(filedir,self.__pdict['runid'])
        vfname='{}ptm_input/dist_velocity_{:04d}.txt'.format(filedir,self.__pdict['runid'])

        # The [x]dfkeys lists give the name of all keys from the appropriate structure that
        # will be written to the corresponding output files and they provide the order
        # in which these will be written.

        # Write the parameter file
        with open(pfname,'w') as pfile:
            # Set the order of parameters to be output
            self.__pfkeys = ['runid','nparticles','ndim','nx','ny','nz','itrace','ifirst','ilast','ntot','dtin','dtout',
                    'istep','iswitch','iphase','nphase','charge','mass','tlo','thi','itraj']
            for key in self.__pfkeys:
                pfile.write('{:<8g} {}\n'.format(self.__pdict[key],key))

        # Write the density file, dist_density_xxxx.dat
        with open(dfname,'w') as dfile:
            # Set the output parameters and their order based on the specified spatial distribution
            if(  self.__ddict['idens']==1):
                dfkeys = ['idens','x0','y0','z0']
            elif(self.__ddict['idens']==2):
                dfkeys = ['idens','xmin','xmax','ymin','ymax','zmin','zmax']
            elif(self.__ddict['idens']==3):
                dfkeys = ['idens','r0','mltmin','mltmax']
            else:
                raise Exception('Density type {:d} is not supported'.format(self.__ddict['idens']))
            for key in dfkeys:
                dfile.write('{:<8g} {}\n'.format(self.__ddict[key],key))

        # Write velocity distribution function file, dist_velocity_xxxx.dat
        with open(vfname,'w') as vfile:
            # Set the output parameters and their order based on the specified velocity distribution
            if(  self.__vdict['idist']==1):
                vfkeys = ['idist','ekev','alpha','phi']
            elif(self.__vdict['idist']==2):
                vfkeys = ['idist','vtperp','vtpara','phi']
            elif(self.__vdict['idist']==3):
                vfkeys = ['idist','nenergy','npitch','phi','emin','emax','pamin','pamax','xsource']
            elif(self.__vdict['idist']==4):
                vfkeys = ['idist','nenergy','npitch','phi','xsource']
            else:
                raise Exception('Distribution type {:d} is not supported'.format(self.__vdict['idist']))
            for key in vfkeys:
                vfile.write('{:<8g} {}\n'.format(self.__vdict[key],key))

        return

    def mlt_to_phi(self,myMlt):

        # Convert from magnetic local time in hours to azimuthal angle in degrees

        if np.isscalar(myMlt):

            if(myMlt < 12.0):
                mlt=myMlt+24
            else:
                mlt=myMlt

        else:

            mlt=myMlt.copy()
            mlt[mlt<12.0]+=24

        res = 15*(mlt-12.0)

        return res

    def create_rungrid(self,times,positions,isSpherical=False):
        """
        -----------
        Description
        -----------
        Given a set of initial positions and start/end times, create the appropriate set of PTM input
        files. A reference text file is also created, which gives the time range and position for each
        runid (the rungrid).

        In the context of this routine, spherical coordinates are [R,MLT,MLAT], where R is in Earth radii, MLT is in
        hours, and MLAT is in degrees.

        ------
        Inputs
        ------
        times(Nt,2)             array or list       First column gives start times, second gives end times
        positions(Nr,3)         array or list       First column gives x-pos, second gives y-pos, third gives z-pos
        isSpherical             logical             Whether coordinates are in spherical (True) or Cartesian (False)
                                                    If true, see note on spherical coordinates above.

        -------
        Outputs
        -------
        None        Creates text files for run inputs and a rungrid file


        ------
        Author
        ------
        Jesse Woodroffe
        jwoodroffe@lanl.gov

        """

        # Particle positions in PTM are Cartesian, but coordinates in the rungrid file are
        # provided in whatever format is sent to this routine.

        if(isSpherical):
            ph=np.deg2rad(np.array([self.mlt_to_phi(mlt) for mlt in positions[:,1]]))
            th=np.deg2rad((90.0-positions[:,2]))
            r=positions[:,0]*np.sin(th)**2
            x=r*np.sin(th)*np.cos(ph)
            y=r*np.sin(th)*np.sin(ph)
            z=r*np.cos(th)
        else:
            x=np.array(positions[:,0])
            y=np.array(positions[:,1])
            z=np.array(positions[:,2])

        tstart = np.array(times[:,0])
        tstop = np.array(times[:,1])

        with open('rungrid.txt','w') as f:
            myid=0
            f.write('{:<8}{:<8}{:<8}{:<8}{:<8}{:<8}\n'.format('RunID','tLo','tHi','x0','y0','z0'))
            for i in range(tstart.size):
                for j in range(x.size):
                    myid+=1
                    f.write('{:<8}{:<8.1f}{:<8.1f}{:<8.3f}{:<8.3f}{:<8.3f}\n'.format(myid,tstart[i],tstop[i],positions[j,0],positions[j,1],positions[j,2]))
                    self.set_parameters(runid=myid,x0=x[j],y0=y[j],z0=z[j],tlo=tstart[i],thi=tstop[i])
                    self.create_input_files()

        return

if __name__ == "__main__":

    """
    This is an example routine that creates a default input object and uses it to make a rungrid
    """

    pic = ptm_input_creator()
    R = 6.6
    MLT = np.arange(24)
    th = (np.pi/180.0)*np.array([rc.mlt_to_phi(lt) for lt in MLT])
    pos = np.zeros([24,3])
    pos[:,0] = R*np.cos(th)
    pos[:,1] = R*np.sin(th)
    tstart = 18000+300*np.arange(48)
    tstop = tstart+7200.0
    times = np.c_[tstart,tstop]

    try:
        os.mkdir('test_rungrid')
        print("Creating test_rungrid subdirectory")
    except:
        print("Existing test_rungrid subdirectory will be used")

    os.chdir('test_rungrid')

    pic.create_rungrid(times,pos)
