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

Last updated: 3/24/2017 Modified create_input_files for greater flexibility using string format methods
rather than c-style string formatters. Removed extraneous global variables

"""

# Define constants for convenience
__ckm = 2.998e5 # Speed of light in km/s

def get_default_inputs(ivelocity=3,idensity=1):
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

  # Change these parameter defaults as appropriate for your application.

  # Default parameters for ptm_parameters file.
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
  pdict['ntot']=2         # Total number of time steps in tgrid file
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
  pdict['itraj']=0        # Write particle trajectories if in flux map mode: 0 = no, 1 = yes

  # Default parameters for velocity file
  vdict={}

  vdict['idist']=ivelocity
  if(ivelocity==1):
    # Ring
    vdict['idist']=1
    vdict['ekev']=100.0     # Particle energy in keV
    vdict['alpha']=90.0     # Particle pitch angle in degrees
    vdict['phi']=180.0      # Particle gyrophase angle in degrees; set negative to seed randomly
  elif(ivelocity==2):
    # Bi-Maxwellian
    vdict['idist']=2
    vdict['vtperp']=0.25*__ckm  # Perpendicular thermal velocity in km/s
    vdict['vtpara']=0.25*__ckm  # Parallel thermal velocity in km/s
    vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
  elif(ivelocity==3):
    # Uniform flux map mode
    vdict['idist']=3
    vdict['nenergy']=34       # Number of energies to use in discrete flux map
    vdict['npitch']=31        # Number of pitch angles to use in discrete flux map
    vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
    vdict['emin']=1.0         # Lowest energy of flux map (keV)
    vdict['emax']=100.0       # Highest energy of flux map (keV)
    vdict['pamin']=0.0        # Lowest pitch angle of flux map (degrees)
    vdict['pamax']=90.0       # Highest pitch angle of flux map (degrees)
    vdict['xsource']=-12.0    # XGSM distance where integration should be terminated (RE)
  elif(ivelocity==4):
    # User-specified flux map mode
    vdict['idist']=3
    vdict['nenergy']=34       # Number of energies to use in discrete flux map
    vdict['npitch']=31        # Number of pitch angles to use in discrete flux map
    vdict['phi']=180.0        # Gyrophase angle in degrees; set negative to seed randomly
    vdict['xsource']=-12.0    # XGSM distance where integration should be terminated (RE)
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

def create_input_files(pf,filedir=''):
  """
  Use a user-provided parameter dictionary to write the ptm input files. The get_default_inputs routine can be used
  to get a template for this dictionary, and individual components can be modified before passing that dictionary
  to this routine.

  * * * * * * * * * *
  ----------
  | INPUTS |
  ----------

  pf    Dictionary    Contains specification of all the parameters for creation of ptm input files:

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

  Last modified 3/24/2017: Simplified routines with new key lists instead of multiple explicit writes. Also
  changed print statements to use string formatting methods.
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
  pfname='{}ptm_parameters_{:04d}.txt'.format(filedir,pf['p']['runid'])
  dfname='{}dist_density_{:04d}.txt'.format(filedir,pf['p']['runid'])
  vfname='{}dist_velocity_{:04d}.txt'.format(filedir,pf['p']['runid'])

  # The xfkeys lists give the name of all keys from the appropriate structure that
  # will be written to the corresponding output files and they provide the order
  # in which these will be written.

  # Write the parameter file
  with open(pfname,'w') as pfile:
      # Set the order of parameters to be output
      pfkeys = ['runid','nparticles','ndim','nx','ny','nz','itrace','ifirst','ilast','ntot','dtin','dtout',
                'istep','iswitch','iphase','nphase','charge','mass','tlo','thi','itraj']
      for key in pfkeys:
        pfile.write('{:<8g} {}\n'.format(pf['p'][key],key))

  # Write the density file, dist_density_xxxx.dat
  with open(dfname,'w') as dfile:
      # Set the output parameters and their order based on the specified spatial distribution
      if(  pf['d']['idens']==1):
        dfkeys = ['idens','x0','y0','z0']
      elif(pf['d']['idens']==2):
        dfkeys = ['idens','xmin','xmax','ymin','ymax','zmin','zmax']
      elif(pf['d']['idens']==3):
        dfkeys = ['idens','r0','mltmin','mltmax']
      else:
        raise Exception('Density type {:d} is not supported'.format(pf['d']['idens']))
      for key in dfkeys:
        dfile.write('{:<8g} {}\n'.format(pf['d'][key],key))

  # Write velocity distribution function file, dist_velocity_xxxx.dat
  with open(vfname,'w') as vfile:
      # Set the output parameters and their order based on the specified velocity distribution
      if(  pf['v']['idist']==1):
        vfkeys = ['idist','ekev','alpha','phi']
      elif(pf['v']['idist']==2):
        vfkeys = ['idist','vtperp','vtpara','phi']
      elif(pf['v']['idist']==3):
        vfkeys = ['idist','nenergy','npitch','phi','emin','emax','pamin','pamax','xsource']
      elif(pf['v']['idist']==4):
        vfkeys = ['idist','nenergy','npitch','phi','xsource']
      else:
        raise Exception('Distribution type {:d} is not supported'.format(pf['v']['idist']))
      for key in vfkeys:
        vfile.write('{:<8g} {}\n'.format(pf['v'][key],key))

  return
