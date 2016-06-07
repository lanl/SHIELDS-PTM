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
  
  None. Three text files are created but the routine returns a null value.
  
  * * * * * * * * * *
  
  See inline commentary for descriptions of individual parameters.
              
  Jesse Woodroffe 
  6/7/2016
  
  """


  # Default parameters for ptm_parameters file
  pdict={}
  pdict['runid']=1
  pdict['nparticles']=1
  pdict['ndim']=3
  pdict['nx']=75
  pdict['ny']=75
  pdict['nz']=75
  pdict['itrace']=-1
  pdict['ifirst']=1
  pdict['ilast']=2
  pdict['dtin']=3600.0
  pdict['dtout']=1.0
  pdict['istep']=1
  pdict['iswitch']=-1
  pdict['iphase']=3
  pdict['nphase']=32
  pdict['charge']=-1.0
  pdict['mass']=1.0
  pdict['tlo']=0.0
  pdict['thi']=3600.0
  
  # Default parameters for velocity file
  vdict={}

  vdict['idist']=ivelocity  
  if(ivelocity==1):
    # Beam
    vdict['idist']=1
    vdict['ekev']=100.0
    vdict['alpha']=90.0
    vdict['phi']=180.0
  elif(ivelocity==2):
    # Maxwellian
    vdict['idist']=2
    vdict['vtperp']=0.25*ckm
    vdict['vtpara']=0.25*ckm
    vdict['phi']=180.0  
  elif(ivelocity==3):
    # Flux map mode
    vdict['idist']=3
    vdict['nenergy']=34
    vdict['npitch']=31
    vdict['phi']=180.0
    vdict['emin']=1.0
    vdict['emax']=100.0
    vdict['pamin']=0.0
    vdict['pamax']=90.0
  else:
    raise Exception('Option ivelocity = '+str(ivelocity)+' is not supported.')
         
  # Default parameters for density file
  ddict={}
  if(idensity==1):
    # Single point
    ddict['idens']=1
    ddict['x0']=-5.0
    ddict['y0']=0.0
    ddict['z0']=0.0 
  elif(idensity==2):
    # Cubic region
    ddict['idensity']=2
    ddict['xmin']=-7.0
    ddict['xmax']=-6.0
    ddict['ymin']=-0.5
    ddict['ymax']= 0.5
    ddict['zmin']=-0.5
    ddict['zmax']=0.5
  elif(idensity==3):
    # Radial distance
    ddict['idens']=3
    ddict['r0']=6.6
    ddict['mltmin']=-6.0
    ddict['mltmax']=6.0    
  
  ptmdict={'p':pdict,'v':vdict,'d':ddict}

  return ptmdict

def create_input_files(pd,filedir=''):
  """
  Use the user-provided parameter dictionary to write the ptm input files.
  
  Jesse Woodroffe
  6/7/2016
  """

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

  dfile=open(filedir+'dist_density_%4.4i.txt' % pd['p']['runid'],'w')
  dfile.write('%-8i \t\t\t\t idens\n' % pd['d']['idens'])
  
  if(pd['d']['idens']==1):
    # Write dist_density file
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

  # Write dist_velocity file
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
    