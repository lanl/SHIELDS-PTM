"""
This file contains routines necessary to create input files for the ptm simulation. See in-line comments
for the description of various parameters.

Usage:

  The simplest case is used to generate a single particle drift trace for dipole field test data:

  config=ptm_input.get_default_inputs()
  ptm_input.create_input_files(config)

  Other input files can be created by modifying the elements of the configuration dictionary, config.
  
Jesse Woodroffe
6/7/2016

"""

from numpy import pi

# Define constants for convenience
ckm = 2.998e5 # Speed of light in km/s
dtor = pi/180.0 # Degrees to radians conversion
rtod = 180.0/pi # Radians to degrees conversion

def get_default_inputs(ivelocity=1,idensity=1):
  """
  Populate the configuration dictionary with appropriate keys and default values. Two input values are used to determine
  what set of parameters are initialized into the dictionary.
  
  * * * * * * * * * *
  ----------
  | INPUTS |
  ----------
  
  This routine has two optional inputs, ivelocity and idensity.   
  
  ivelocity: Integer [1 2 3]
  
  1 = Monoenergetic beam (default); sets the following parameters:
    E0    Beam energy (keV)
    A0    Beam pitch angle (degrees)
    PHI0  Particle gyrophase angle (degrees), set negative to seed randomly
  2 = Maxwellian distribution; sets the following parameters:
    VTPERP  Perpendicular thermal velocity (km/s)
    VTPARA  Parallel thermal velocity (km/s)
    PHI0    Particle gyrophase angle (degrees), set negative to seed randomly
  3 = Flux map mode; sets the following parameters:
    NENERGIES     Number of energies to use
    NPITCHANGLES  Number of pitch angles to use
    PHI0          Particle gyrophase angle (degrees), set negative to seed randomly
    EKEVMIN       Lowest energy of flux map (keV)
    EKEVMAX       Highest energy of flux map (keV)
    PAMIN         Minimum pitch angle of flux map (degrees)
    PAMAX         Maximum pitch angle of flux map (degrees)
            
  idensity: Integer [1 2 3]

  1 = All particles at same position (default); sets the following parameters:
    X0  Initial x-position (RE)
    Y0  Initial y-position (RE)
    Z0  Initial z-position (RE)
  2 = Randomly seed throughout rectangular region; sets the following parameters:
    XMIN  Minimum x-limit of region (RE)
    XMAX  Maximum x-limit of region (RE)
    YMIN  Minimum y-limit of region (RE)
    YMAX  Maximum y-limit of region (RE)
    ZMIN  Minimum z-limit of region (RE)
    ZMAX  Maximum z-limit of region (RE)
  3 = Randomly seed at fixed geocentric distance; sets the following parameters:
    R0      Radial distance of seeding (RE)
    MLTMIN  Minimum magnetic local time of seeding (HOURS)
    MLTMAX  Maximum magnetic local time of seeding (HOURS)
  
  -----------
  | OUTPUTS |
  -----------
  
  ptmdict: Dictionary
  
  ptmdict has numerous keys that correspond to different parameters used to configure the ptm simulation. See
  the code itself for descriptions of individual keys.
  
  * * * * * * * * * *
  
  See inline commentary for descriptions of individual parameters.
              
  Jesse Woodroffe 
  6/7/2016
  
  """

  # Default parameters for ptm_parameters file
  pdict={}
  pdict['runid']=1
  pdict['nparticles']=1   # Number of particles to use in simulation
  pdict['ndim']=3         # Number of spatial dimensions of input data (2 or 3)
  pdict['nx']=75          # Number of points in X grid
  pdict['ny']=75          # Number of points in Y grid
  pdict['nz']=75          # Number of points in Z grid
  pdict['itrace']=-1      # 1 = forward in time, -1 = backwards in time
  pdict['ifirst']=1       # Timestep of first fields snapshot: e.g., start at bx3d_0005.bin if ifirst=5
  pdict['ilast']=2        # Timestep of last fields snapshot: e.g., finish at bx3d_0009.bin if ilast=9
  pdict['dtin']=3600.0    # Time between snapshots (seconds). If dtin<0, read times from tgrid.bin
  pdict['dtout']=1.0      # Time between data outputs (seconds)
  pdict['istep']=1        # Integrator to use: 1=RK4, 2=Adaptive RKSuite
  pdict['iswitch']=-1     # Particle equations to solve: -1=guiding center, 0=dynamic switching, 1=full orbit
  pdict['iphase']=3       # Gyrophase switching method: 1=random, 2=brute force search, 3=gradient method
  pdict['nphase']=32      # Number of points for brute force search (minimize difference in B)
  pdict['charge']=-1.0    # Particle charge in multiples of fundamental
  pdict['mass']=1.0       # Particle mass in multiples of the electron mass
  pdict['tlo']=0.0        # Lowest time of particle trace
  pdict['thi']=3600.0     # Highest time of particle trace
  
  # Default parameters for velocity file
  vdict={}

  vdict['idist']=ivelocity  
  if(ivelocity==1):
    # Beam
    vdict['idist']=1
    vdict['ekev']=100.0     # Particle energy in keV
    vdict['alpha']=90.0     # Particle pitch angle in degrees
    vdict['phi']=180.0      # Particle gyrophase angle in degrees; set negative to seed randomly
  elif(ivelocity==2):
    # Maxwellian
    vdict['idist']=2
    vdict['vtperp']=0.25*ckm  # Perpendicular thermal velocity in km/s
    vdict['vtpara']=0.25*ckm  # Parallel thermal velocity in km/s
    vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
  elif(ivelocity==3):
    # Flux map mode
    vdict['idist']=3
    vdict['nenergy']=34       # Number of energies to use in discrete flux map
    vdict['npitch']=31        # Number of pitch angles to use in discrete flux map
    vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
    vdict['emin']=1.0         # Lowest energy of flux map (keV)
    vdict['emax']=100.0       # Highest energy of flux map (keV)
    vdict['pamin']=0.0        # Lowest pitch angle of flux map (degrees)
    vdict['pamax']=90.0       # Highest pitch angle of flux map (degrees)
  else:
    raise Exception('Option ivelocity = '+str(ivelocity)+' is not supported.')
         
  # Default parameters for density file
  ddict={}
  if(idensity==1):
    # Single point
    ddict['idens']=1
    ddict['x0']=-5.0      # Initial X-position (RE)
    ddict['y0']=0.0       # Initial Y-position (RE)
    ddict['z0']=0.0       # Initial Z-position (RE)
  elif(idensity==2):
    # Cubic region
    ddict['idensity']=2
    ddict['xmin']=-7.0    # Minimum X-boundary of seeding region (RE)
    ddict['xmax']=-6.0    # Maximum X-boundary of seeding region (RE)
    ddict['ymin']=-0.5    # Minimum Y-boundary of seeding region (RE)
    ddict['ymax']= 0.5    # Maximum Y-boundary of seeding region (RE)
    ddict['zmin']=-0.5    # Minimum Z-boundary of seeding region (RE)
    ddict['zmax']=0.5     # Maximum Z-boundary of seeding region (RE)
  elif(idensity==3):
    # Radial distance
    ddict['idens']=3
    ddict['r0']=6.6       # Radial distance for seeding (RE)
    ddict['mltmin']=-6.0  # Minimum local time for seeding region (HOURS)
    ddict['mltmax']=6.0   # Maximum local time for seeding region (HOURS)
  
  ptmdict={'p':pdict,'v':vdict,'d':ddict}

  return ptmdict

def create_input_files(pd,filedir=''):
  """
  Use a user-provided parameter dictionary to write the ptm input files. The get_default_inputs routine can be used
  to get a template for this dictionary, and individual components can be modified before passing that dictionary
  to this routine.
  
  * * * * * * * * * *
  ----------
  | INPUTS |
  ----------  
  
  pd    Dictionary    Contains specification of all the parameters for creation of ptm input files:
    
  1. ptm_parameters_xxxx.dat
  2. dist_velocity_xxxx.dat
  3. dist_density_xxxx.dat
  
  filedir   String  Optional    Directory (relative or absolute) where output files should be written.
  
  See description of individual parameters in <get_default_inputs>

  -----------
  | OUTPUTS |
  -----------  
  
   None. Three text files are created but the routine returns a null value.

  * * * * * * * * * * 
  
  Jesse Woodroffe
  6/7/2016
  """

  # Handle a couple of special cases
  if(len(filedir)>0):
    if(fileir=='.'):
      # User specified present directory using dot; strip the dot 
      fileir=''
    elif(searchDir[-1]!='/'): 
      # User forgot to add a slash at the end of the directory name; add a slash
      filedir+='/'

  # Write the parameters file, ptm_parameters_xxxx.dat
  
  pfile=open(filedir+'ptm_parameters_%4.4i.txt' % pd['p']['runid'],'w')
  pfile.write('%-8i \t\t\t\t RunID this needs to match the 4 digit tag (without leading zeros)\n' % pd['p']['runid'])
  pfile.write('%-8i \t\t\t\t Number of particles (value is overridden for dist_velocity::idist==3)\n' % pd['p']['nparticles'])
  pfile.write('%-8i \t\t\t\t Spatial dimensions of the fields, 2 or 3: if 2, NZ is automatically set to 2\n' % pd['p']['ndim'])
  pfile.write('%-8i \t\t\t\t Number of cells in x-direction\n' % pd['p']['nx'])
  pfile.write('%-8i \t\t\t\t Number of cells in y-direction\n' % pd['p']['ny'])
  pfile.write('%-8i \t\t\t\t Number of cells in z-direction\n' % pd['p']['nz'])
  pfile.write('%-8i \t\t\t\t Direction of trace, -1 = backwards from end, +1 = forwards from start\n' % pd['p']['itrace'])
  pfile.write('%-8i \t\t\t\t Index of first data file\n' % pd['p']['ifirst'])
  pfile.write('%-8i \t\t\t\t Index of last data file\n' % pd['p']['ilast'])
  pfile.write('%-8.1f \t\t\t\t Cadence of input files (s), if negative get values from tgrid\n' %pd['p']['dtin'])
  pfile.write('%-8.1f \t\t\t\t Cadence of output files (s)\n' % pd['p']['dtout'])
  pfile.write('%-8i \t\t\t\t Time integrator, 1= RK4, 2= RKSuite\n' % pd['p']['istep'])
  pfile.write('%-8i \t\t\t\t Equations to solve: -1 = drift, 0 = switch, 1 = orbit\n' % pd['p']['iswitch'])
  pfile.write('%-8i \t\t\t\t Gyrophase switching method 1 = random, 2 = brute force, 3 = gradient\n' % pd['p']['iphase'])
  pfile.write('%-8i \t\t\t\t Number of brute force points to search (angular resolution = 360/nphase) \n' % pd['p']['nphase'])
  pfile.write('%-8.1f \t\t\t\t charge in multiples of fundamental\n' % pd['p']['charge'])
  pfile.write('%-8.1f \t\t\t\t mass in multiples of electron mass\n' % pd['p']['mass'])
  pfile.write('%-8.1f \t\t\t\t lower limit of time integration (s)\n' % pd['p']['tlo'])
  pfile.write('%-8.1f \t\t\t\t upper limit of time integration (s)\n' % pd['p']['thi'])
  pfile.close()

  # Write the density file, dist_density_xxxx.dat

  dfile=open(filedir+'dist_density_%4.4i.txt' % pd['p']['runid'],'w')
  dfile.write('%-8i \t\t\t\t idens\n' % pd['d']['idens'])
  
  if(pd['d']['idens']==1):
    dfile.write('%-8.2f \t\t\t\t Initial X position (Re)\n' % pd['d']['x0'])
    dfile.write('%-8.2f \t\t\t\t Initial Y position (Re)\n' % pd['d']['y0'])
    dfile.write('%-8.2f \t\t\t\t Initial Z position (Re)\n' % pd['d']['z0'])
  elif(pd['d']['idens']==2):
    dfile.write('%-8.2f \t\t\t\t Minimum X Position (Re)\n' % pd['d']['xmin'])
    dfile.write('%-8.2f \t\t\t\t Maximum X Position (Re)\n' % pd['d']['xmax'])
    dfile.write('%-8.2f \t\t\t\t Minimum Y Position (Re)\n' % pd['d']['ymin'])
    dfile.write('%-8.2f \t\t\t\t Maximum Y Position (Re)\n' % pd['d']['ymax'])
    dfile.write('%-8.2f \t\t\t\t Minimum Z Position (Re)\n' % pd['d']['zmin'])
    dfile.write('%-8.2f \t\t\t\t Maximum Z Position (Re)\n' % pd['d']['zmax'])
  elif(pd['d']['idens']==3):
    dfile.write('%-8.2f \t\t\t\t Radial Distance (Re)\n' % pd['d']['r0'])
    dfile.write('%-8.2f \t\t\t\t Minimum MLT (Re)\n' % pd['d']['mltmin'])
    dfile.write('%-8.2f \t\t\t\t Maximum MLT (Re)\n' % pd['d']['mltmax'])  
  else:
    raise Exception('Density type '+str(pd['d']['idens'])+' is not supported')
  
  dfile.close()

  # Write velocity distribution function file, dist_velocity_xxxx.dat

  vfile=open(filedir+'dist_velocity_%4.4i.txt' % pd['p']['runid'],'w')
  vfile.write('%-8i \t\t\t\t idist\n' % pd['v']['idist'])
  if(pd['v']['idist']==1):
    vfile.write('%-8i \t\t\t\t E(keV)\n' % pd['v']['ekev'])
    vfile.write('%-8i \t\t\t\t Pitch Angle (deg) \n' % pd['v']['alpha'])
    vfile.write('%-8.2f \t\t\t\t Phase Angle\n' % pd['v']['phi'])    
  elif(pd['v']['idist']==2):
    vfile.write('%-8i \t\t\t\t Maxwellian vtPerp\n' % pd['v']['vtperp'])
    vfile.write('%-8i \t\t\t\t Maxwellian vtPara\n' % pd['v']['vtpara'])
    vfile.write('%-8.2f \t\t\t\t Phase Angle\n' % pd['v']['phi'])        
  elif(pd['v']['idist']==3):
    vfile.write('%-8i \t\t\t\t nEnergies\n' % pd['v']['nenergy'])
    vfile.write('%-8i \t\t\t\t nPitchAngles\n' % pd['v']['npitch'])
    vfile.write('%-8.2f \t\t\t\t Phase Angle\n' % pd['v']['phi']) 
    vfile.write('%-8.2f \t\t\t\t Emin(keV)\n' % pd['v']['emin'])
    vfile.write('%-8.2f \t\t\t\t Emax(keV)\n' % pd['v']['emax'])
    vfile.write('%-8.2f \t\t\t\t Alphamin(deg)\n' % pd['v']['amin'])
    vfile.write('%-8.2f \t\t\t\t Alphamax(deg)\n' % pd['v']['amax'])  
  else:
    raise Exception('Distribution type '+str(pd['v'][idist])+' is not supported')

  vfile.close()
  
  return
    