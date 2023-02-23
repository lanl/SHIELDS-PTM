# SHIELDS-PTM: Particle Tracing Model

## Installing SHIELDS-PTM
SHIELDS-PTM is developed under version control at https://github.com/lanl/SHIELDS-PTM
Tagged releases are also archived at Zenodo. The latest release is always accessible via
https://doi.org/10.5281/zenodo.4891973

Once the desired version of the code is obtained, simply follow the instructions in ```README.md```.
The python module can also be installed using the standard ```python setup.py install``` if desired.
Test status for the SHIELDS-PTM code is displayed at the github repository, but to ensure that the
code runs locally, we provide a simple test case. After compilation, simply run ```./ptm 1```
To use the ```ptm_python``` module, simply ```import ptm_python``` in your interpreter or script.

## Using SHIELDS-PTM
### Preparing PTM input files

Python functions and classes used to prepare PTM input files are in the ```ptm_python``` module.
Some more specialized, or targeted, scripts for specific applications are in the ```scripts```
directory.

For a run, a run directory should be created. This can have any name and should contain directories
named ```ptm_data```, ```ptm_input```, and ```ptm_out```. It should also contain a soft-link to the
```ptm``` executable.

The run directory should be set up with the following structure:
```
.
+-- rundir
    |   `-- ptm
    +-- ptm_input
    |   |-- ptm_parameters_rrrr.txt
    |   |-- dist_density_rrrr.txt
    |   `-- dist_velocity_rrrr.txt
    +-- ptm_data
    |   |-- ptm_fields_nnnn.txt
    |   |-- tgrid.dat
    |   |-- energies_rrrr.bin
    |   `-- pitchangles_rrrr.bin
    `-- ptm_output
```
The ```rrrr``` is a 4-digit, zero-padded _run_ number, in case the simulation is partitioned into multiple runs.
The ```nnnn``` is a 4-digit, zero-padded _timestep_ for the input fields file. Note that for a static field, two identical timesteps should be supplied.

- If using SWMF/BATSRUS electric and magnetic fields, use ```ptm_python/ptm_preprocessing.py``` as a script to generate the fields and time grid files, and note that the preferred format is the 3D "IDL" files. The preprocessing script can also use the ```3D MHD tec``` output from BATSRUS.
- If using semi-empirical magnetic field models, use the ```ptm_lanlgeomag``` command-line program in the ```tools``` directory. This requires that the LANLGeoMag library is installed. LANLGeoMag is available from https://github.com/drsteve/LANLGeoMag

Full details of the file specifications are given in the File Formats section of this document.

#### ptm_preprocessing.py

2. calls "ptm_read.py" to read input ascii files ".mhd" (at all n epochs).
   The Tecplot format files must use a specific output format (currently not documented, probably `3d MHD tec`)
   A nominal column layout is given here:
    x   y   z   rho    ux   uy   uz    bx  by  bz   P    jx     jy     jz
    RE  RE  RE  Mp/cc  km/s km/s km/s  nT  nT  nT   nPa  muA/m2 muA/m2 muA/m2
   The IDL format files use the flexible reader from SpacePy and hence any supported format
   can be used.
3. calls function "gauss_interp_EB" of module "ptm_interpolate.py"
   a. to interpolate B from the SWMF mesh (of .mhd) to the uniform XYZ grid created above (1b)
   b. to calculate electric field E=-uxB and write files "e(x,y,z)3d_????.bin"
1  c. write files "b(x,y,z)3d_????.bin" at all n input epochs

Module ```ptm_input.py``` will generate run configuration files in directory ptm_input/
 a. file "ptm_parameters_????.txt" (with ????=4-digit zero-padded runid above) containing global pars:
      number of particles, number of nx,ny,nz cells, number of timesteps
      cadence in output trajectory file ptm_output/ptm_????.dat
      cadence of input E&B files in ptm_data (note that the 4-digit number on these files
      indicates the timestep, and these files are used by ALL runs)
      indices ifirst and ilast of first/last epochs of E&B files to be read by PTM, etc.
    IMPORTANT:  Choose these global pars and make sure that the spatial grid has same number of nodes 
                (nx)(ny)(nz) as outputted by module ```ptm_preprocessing.py``` in ```ptm_parameters``` file (see 1b above).
 b. file "dist_density_????.txt" with parameters for the spatial distribution:
    first line:
      idens=1 -> single-location in RE
      idens=2 -> cube corners in RE
      idens=3 -> fixed radial distance in RE in equatorial plane, MLT range in hours
 c. file "dist_velocity_????.txt" with parameters for the energy/PA distributions:
    first line:
      idist=1 -> single E/keV, single PA/deg, single or random phase angle
      idist=2 -> vtperp and vtpara in km/s, single or random phase angle
      idist=3 -> uniform flux-map mode: # of Energies, # of PitchAngles, gyropase phi, low/high E/PA
      idist=4 -> user-specified flux map mode: as above (first three)

 IMPORTANT: For each of the 3 switches above (runid, idens, idist), set all parameters (after -> above).               
  Alternatively, you can just edit existing files in dir ptm_input/ and remember that PTM will read the 
  three input files that have the N=runid specified when PTM is run.

--------------

    NOTE: The two scripts, ptm_preprocessing.py and ptm_input.py, for preparing the PTM input files have in common
         only one thing: the number of XYZ grid nodes nx*ny*nz, and you must use the same nx,ny,nz in both codes 
         or, else, PTM will crash, get to "X/Y/Z out of bounds" or churn out lots of NaN's, without telling you why.
         See A2.1.b and A3.a/NOTE above.

-----------------------------------------------

### Running PTM

 I. INPUT 
    The SHIELDS-PTM simulation is configured using following files:

 1a. General parameters are in file "ptm_input/ptm_parameters_????.txt"
   called by subroutine "read_ptm_parameters" in module FILEIO
   sets global parameters for simulation (see above)

 1b. Particle spatial distribution is in file "ptm_input/dist_density_????.txt"
   called by subroutine "particle_initialize" in module PARTICLES
   (particle location goes into object myParticle%x)
   cases on first line:
     idens=1 -> single-location
     idens=2 -> random seeding in a cube
     idens=3 -> random MLT and fixed radial distance in equatorial plane
   there is no user-specified spatial distribution

 1c. Particle energy/PA distribution is in file "ptm_input/dist_velocity_????.txt"
   called by subroutine "particle_initialize" in module PARTICLES
   (to initialize vpara and vperp distributions of object myParticle%v in Cartesian XYZ coordinates,
      for use in orbit equations (iswitch=1 in file ptm_parameters))
   cases on first line:
      idist=1 -> monoenergetic particles, single PA
      idist=2 -> Maxwellian distribution
      idist=3 -> uniform distribution in E and PA (uniform flux map mode)
      idist=4 -> user-specified flux map mode (Es and PAs in input files energies.bin and pitchangles.bin)

 1d. Time/Space grids and EM field input data (at least two epochs) are in dir "ptm_data/"
    files read with interface "read_array" in module GLOBAL (data stored in "C-style")

    1d1. Files "ptm_data/x,y,zgrid.bin" contain linear spatial grid read by subroutine "fields_initialize" in module FIELDS
       In file "ptm_parameters", dtIn is the cadence of E&B input files.
       If dtIn>0, then time-grid (for E&B) is constructed linearly between Tmin=0 and Tmax=dtIn*nt, 
         with nt = ilast-ifirst (from ptm_parameters) = number of E&B input files
       If dtIn=<0, then time-grid (for E&B) Tmin-Tmax is read from file "ptm_data/tgrid.bin" 

    1d2. Files "ptm_data/(b,e)(x,y,z)3d_????.bin" contain E&B at epochs from N=ifirst to N=ilast (see ptm_parameters)
       They are read by subroutine "fields_initialize" in module FIELDS
       Subroutine "get_fields" in module FIELDS calculates B and E
         at any point in the simulation cube through interpolation (calls subroutine tricubic_interpolate in module INTERPOLATION)
         at each timestep though linear interpolation

<!-- pagebreak -->
## Numerical Approach

 1. Guiding-center (GC) drift mode vs Full-orbit (FO) integration

 Particle time-advancing is done in subroutine "stepper_push" in module STEPPER, which integrates the particle equation
 of motion either in the GC representation (logical "myParticle%drift") or in the FO mode.
 Subroutine "particle_initialize" initializes the particle trajectory integration in FO mode.

 IF iswitch=-1 (in "ptm_parameters"), then only GC mode is used. Validity of GC mode (local B-field uniformity - below) is not verified.
 IF iswitch=1,                        then only FO mode is used. For non-rel electrons, time-advancing may be much less than requested.
 IF iswitch=0, then the integration method is switched between GC and FO depending on the local magnetic uniformity:
  a. If already in GC mode and the length-scale/curvature radius of magnetic field B/gradB (~RE ?) is smaller than
     Rgyr/kappa_orbit=100*Rgyr, then integration is switched to FO mode.
  b. If already in FO mode and the length-scale/curvature radius of magnetic field B/gradB is larger than
     Rgyr/kappa_drift=125*Rgyr, then integration is switched to GC mode.

  Note: gyration radius Rgyr=4.4km*(250nT/B)*sin(PA)*sqrt(E/100keV) for non-rel els, a factor sqrt(mp/me)=43 larger for pr -> 180 km
                        Rgyr=13.3km*(250nT/B)*sin(PA)*(E/1MeV)      for rel     els,
     thus electrons of E < 5 MeV and protons of E < 15 keV satisfy the magnetic field uniformity criterion: 100*Rgyr < B/gradB ~ RE,
     thus, for energies of interest, electrons stay in GC mode, protons stay in FO mode.

  The GC/FO switch is done with a choice of gyrophase, based on the "iphase" parameter in file "ptm_parameters":
   iphase=1 -> random gyrophase, iphase=2 -> simple minimization of delta B, iphase=3 -> Pfefferle's method (to minimize delta B).

 The GC mode allows larger timesteps than FO mode, thus GC mode
  a. leads to shorter simulations
  b. conserves better particle energy


 2. Adaptive timestep

 Object myParticle is initialized in subroutine "particle_initialize" of module PARTICLES.
 Subroutine "stepper_push" advances a particle with smaller timesteps adding up to the output-file cadence dtOut,
 by updating quantities (position, velocity) of object myParticle in either GC or FO representations.
 The integration timestep of "stepper_push" is set initially in subroutine "particle_initialize".
 Subroutine "stepper_push" re-calculates timestep duration dt after each time advance
   a. IF in FO mode: dt is a small fraction (GLOBAL parameter epsilon_orbit=0.05) of the orbital period/(2*pi)
   b. IF in GC mode: dt is a small fraction (GLOBAL parameter epsilon_drift=0.1) of the timescale at which all accelerations
         involved (electric, gradient, curvature, mirror) would yield the current particle speed (the "Kress-Hudson heuristic")

 The GC mode allows larger timesteps than FO mode, thus GC mode
  a. leads to shorter simulations
  b. conserves better particle energy
  From 1. NOTE above: for energies of interest (E_el < 1 MeV and E_pr > 100 keV), electrons are in GC mode, protons in FO mode.


 3. PTM and STEPPER time-advancing

 The total number of timesteps is set in modules PTM and STEPPER.
 PTM calls subroutine "stepper_push" to advance the particle by the output-file cadence dtOut (set in file ptm_parameters).
 After each dtOut, particle position and velocity (converted to lab-frame from GC or FO) are stored for output.
 PTM outputing stops when the entire simulation time Thi (set in ptm_parameters) has been reached or when its time-advancing
  loop is stopped by "stepper_push" (for reasons below).
 In PTM, the number of time-steps and output calls is nwrite = (Thi-Tlo)/dtOut.

 In STEPPER, time-advancing is stopped if the drift or orbit integrator encounter a problem, if the particle hits the ionosphere
 (r=1 R_E), or if the particle has reached the source region boundary (when in flux map mode).
  Because time-advancing is prescribed to be done a fixed number (nstep_max) of times, based on the initial timestep duration
 (dt_base), the adaptive timestep may lead to a time-advancing less than that requested by PTM (dtOut).
  This is likely to happen for electrons in FO mode, where the timestep is 1/100 of the gyration period, because of their small
 Tgyr, but should not be a problem for protons because Tgyr ~ mass (mass is a parameter in file "ptm_paramaters") or for any
 particle in GC mode.


 4. Timestepping

 4A. In GC mode
 Timestep size dT is set by particle velocity vel and all accelerations involved acc.
 Consequently, for a fixed advancing time Thi, run duration Trun depends on particle PA, energy E, and mass m
 (besides E and B fields): Trun ~1/dT ~vel/acc. If acc~force/m, then a naive expectation would be that
 Trun ~ 1/sqrt(mE) (after vel~sqrt(E/m)), but not accelerations satisfy acc~1/m.

 Empirically, run-time dependence on
  a. PA     is = Trun(PA=90)*(10-PA/10) - longer runs for smaller PA
  b. Energy is ~ E^0.6     at PA=90     - longer runs for higher energy (exponent decreases with PA)
  c. mass   is ~ m^{-0.4}               - longer runs for smaller mass
      (but this is not so relevant because only electrons should be in GC mode,
       while next lightest particle - protons - should be in FO mode).
 So, Trun ~ sqrt(E/m) --> expectation Trun~sqrt{E} is wrong
                      --> expectation Trun~1/sqrt{m} is right
       (despite that accelerations in GC frame are mass-independent !)

 Example: for B=250 nT, a (timestep propto) parameter epsilon_drift = 0.1 (set in module GLOBAL) requires
  a. ~50 steps to advance an Electron for dtOut=1s, leading to Trun=10s for a Thi=1h simulation
  b. few steps to advance a  Proton   for dtOut=1s, leading to Trun=.5s for a Thi=1h simulation

 Note: Electron trajectories should be calculated in GC mode because electrons of energy E < 10 MeV satisfy
       the GC mode's magnetic uniformity criterion (see above) and GC mode timesteps are larger than for FO,
       so runs are shorter.


 4B. In FO mode

 Testing the full orbital motion using a 100keV electron at L=5 (pitch angle ~ 90deg), one timestep takes about 50 microsecs on a relatively modern test machine.

 Timestep size dT is set by particle gyration period T_g ~ gamma*m/B, thus run duration,
  a. depends on particle mass Trun ~ 1/m  (relevant for ions, which do not satisfy the GC mode criterion),
  b. is independent on Energy if particle is non-relativistic (otherwise T_g ~ E and Trun ~ 1/E),
  c. is dependent on PA through T_g ~ 1/B, with B along the field line depending on (equatorial) PA.
     Empirically, Trun(PA) = Trun(PA=90)*(0.6+35/PA) - longer runs for smaller PA

 Example: for B=250 nT, a (timestep propto) parameter epsilon_orbit = 0.05 (set in module GLOBAL) requires
  a. ~2/3 million steps to advance an Electron for dtOut=1s, leading to Trun=20h  for Thi=1h
  b. ~300         steps to advance a  Proton   for dtOut=1s, leading to Trun=1min for Thi=1h

 Note: Proton trajectories should be calculated in FO mode because protons of energy E > 0.1 MeV do not
       satisfy GC mode's magnetic uniformity criterion.



 5. Integration of ODE for particle motion

 In subroutine ```push```, time-integration can be done with either
  - a fixed timestep, 4th-order Runge-Kutta integrator: subroutine ```rk4``` in module ```pusher```, if istep=1 in file "ptm_parameters" or
  - an adaptive timestep Runge-Kutta integrator with error control: subroutine ```range_integrator``` in module ```rksuite```, if istep=2 in file "ptm_parameters"
 without a dynamical switching between these integrators.
 Note: The RKSuite adaptive integrator is much more computationally-expensive than the fixed-step RK4 method.


 6. Particle periodic motions

 "ptm_????.dat" trajectory file has cadence dtOut=1 sec (set in file "ptm_parameters_????.txt").
 a. Gyration period is Tgyr = 0.14ms*(250nT/B) for electrons, a factor mp/me=1837 larger for protons -> 0.26s
     (independent of particle energy E for non-relativistic particles)
    Conclusion: particle orbit is (much) under-sampled in output file "ptm_????.dat"
 b. Bounce period is Tb = (0.6-0.9)s x sqrt(100keV/E)(R/5RE) for non-relativistic electrons, a factor sqrt(mp/me)=43 larger for protons
     coefficient is weakly dependent on pitch angle because mirror point goes to a larger zmax with decreasing PA but vpara
      increases with decr PA, so that zmax/cos(PA) is "stiff", from PTM: coeff=0.58 for PA=80, coef=0.93 for PA=10
    Conclusion: bounce motion is under-sampled for els, well-sampled for prot in output file "ptm_N.dat"
 c. Particle drift is dominated by radial dependence of B (ExB drift has a velocity 10x smaller).
    Drift period Td = 1.8h(B/250nT)(100keV/E)(R/5RE)^2 (any particle) has a very weak dependence on PA
     (note drift velocity ~ vperp^2/B = const along field line, due to 1st invariant -> spatial drift constant during bounce
      factor 1: smaller PA implies smaller vperp in equatorial plane and smaller spatial drift velocity
      factor 2: smaller PA implies access to higher latitudes, where the angular drift is larger for constant drift velocity,
      1 & 2 compensate each other).
  Conclusion: drift motion is well-sampled in output file "ptm_????.dat"


## Output Files

 The PTM simulation outputs in two different types of output files:
 "ptm_????.dat" if idist=1,2 (monoenergetic particles/Maxwellian distribution)
 "map_????.dat" if idist=3,4 (uniform E&PA distributions/user-specified flux map mode)
 (idist from "dist_velocity" file)
 The ptm trajectory files can be enabled for fluxmap mode.

 The "ptm" output file contains trajectories for each particle in the simulation,
 Time    X-Pos   Y-Pos   Z-Pos   gamma*Vperp   gamma*Vpara   Energy   PitchAngle
 note: velocities are v_ptm = gamma*v_real, with gamma=(Ek/mc2+1)=1/sqrt(1-vreal^2)=sqrt(1+vptm^2) !!!

 When SHIELDS-PTM is run in parallel, particles are not output in any particular order. 
 Instead, a particle's full time history is output as soon as it is available (i.e. as soon as the particle has 
 finished integration). Because of the inherent disorder of the "ptm" file, it is not convenient to use this file 
 for Liouville tracing estimates of the distribution, which is one of the primary reasons that SHIELDS-PTM was 
 developed in the first place.
 
 The "map" file is designed for this purpose, and lists (for each particle):
  Time  X-Pos   Y-Pos   Z-Pos (at end)
  Energy  PitchAngle          (at start)
  Energy                      (at end)
  Velocity (Cartesian)        (at start)
 (run-end is initial-time if itrace = -1 = trace backwards from start, or final-time if itrace = 1 = trace forward from start)
 Trajectories are, obviously, not recorded.

 UNITS:
 Location/distance in R_E
 velocity in km/s
 energy in keV
 B mag in nT
 E in nT*km/s= muV/m



----------------------------------------

 C. Processing of PTM output

 Module ptm_postprocessing.py processes PTM output fluxmap files
   to produce differential and omni-directional fluxes on the grid of Energies and PitchAngles
   created by PTM with options idist=3 (uniform grid) or idist=4 (user-specified grid) (see above).

 Method "ptm_postprocessing.process_run" calls
  1. "parse_map_file" of module PTM_TOOLS.py" to read the "map" file
       and create a raw-flux dictionary fluxmap, which is further used by
  2. "calculate_flux" that yields differential flux for an assumed initial particle kappa-distribution 
      (choose index kappa and charcateristic energy e_char), further used by
  3. "calculate_omni-directional_flux" that calculates omni-flux,
  and outputs a dictionary containing: 
  1a. Energy x PA array = number of particles
  1b. final_energy (E,PA) - particle energy at simulation end (see itrace)
  1c. final_position-x,y,z (E,PA) - particle location at end
  2a. Energy grid 
  2b. PA grid
  2c. differential flux (E,PA)
  2d. omni-flux (E)

 UNITS:
  flux in 1/(keV cm^2 s sr)
  omni-flux in 1/(keV cm^2 s)

## File Formats

### Input files: PTM Fields
As noted previusly, several SWMF/BATSRUS output formats are supported by the SHIELDS-PTM preprocessing.
The 3D "IDL" format files use the flexible SpacePy reader so most variable sets can be used (SI unit output is required).
The Tecplot format files from BATSRUS must use a specific output format (```3d MHD tec```). A nominal column layout is given here:
```
x   y   z   rho    ux   uy   uz    bx  by  bz   P    jx     jy     jz
RE  RE  RE  Mp/cc  km/s km/s km/s  nT  nT  nT   nPa  muA/m2 muA/m2 muA/m2
```


#### Notes on format
The file is ASCII and will have a filename following the
ptm_fields_XXXX.dat format, where XXXX is a zero-padded number
indicating the timestep.

#### Header
The first line of the header gives the (integer) number of grid
points in each dimension, e.g., for a 75x75x75 grid:
```
75 75 75
```

The next 3 lines give the coordinates of the grid points. First
the X-coordinates are given, followed by a newline. Then the Y
and Z coordinates. For a regular cube in (X,Y,Z) centered on the
origin, these three lines will be identical. For a 75x75x75 grid
centered at (0, 0, 0), extending from -6 to 6 in each dimension,
the lines will be as follows:
```
-6.00000e+00 -5.83784e+00 -5.67568e+00 -5.51351e+00 -5.35135e+00 -5.18919e+00 -5.02703e+00 -4.86486e+00 -4.70270e+00 -4.54054e+00 -4.37838e+00 -4.21622e+00 -4.05405e+00 -3.89189e+00 -3.72973e+00 -3.56757e+00 -3.40541e+00 -3.24324e+00 -3.08108e+00 -2.91892e+00 -2.75676e+00 -2.59459e+00 -2.43243e+00 -2.27027e+00 -2.10811e+00 -1.94595e+00 -1.78378e+00 -1.62162e+00 -1.45946e+00 -1.29730e+00 -1.13514e+00 -9.72973e-01 -8.10811e-01 -6.48649e-01 -4.86486e-01 -3.24324e-01 -1.62162e-01  0.00000e+00  1.62162e-01  3.24324e-01  4.86486e-01  6.48649e-01  8.10811e-01  9.72973e-01  1.13514e+00  1.29730e+00  1.45946e+00  1.62162e+00  1.78378e+00  1.94595e+00  2.10811e+00  2.27027e+00  2.43243e+00  2.59459e+00  2.75676e+00  2.91892e+00  3.08108e+00  3.24324e+00  3.40541e+00  3.56757e+00  3.72973e+00  3.89189e+00  4.05405e+00  4.21622e+00  4.37838e+00  4.54054e+00  4.70270e+00  4.86486e+00  5.02703e+00  5.18919e+00  5.35135e+00  5.51351e+00  5.67568e+00  5.83784e+00  6.00000e+00
```

All subsequent rows are data in the following format:
```
i j k Bx[i,j,k] By[i,j,k] Bz[i,j,k] Ex[i,j,k] Ey[i,j,k] Ez[i,j,k]
```
where i, j and k are the index into the grid coordinates.
The indexing starts from 1 (not zero), thus for a 75x75x75 grid
the rows will start with indices
```
   1    1    1 ...
   1    1    2 ...
   1    1    3 ...
   .
   .
   .
   1    1   75 ...
   1    2    1 ...
   .
   .
   .
   75  75  74 ...
   75  75  75 ...
```
All coordinate indices are integer, and all cordinate and data values
are floating point. The usual format is `[-]1.23456e+00`, that is,
six digits plus a two digit exponent. The sign is optionally included
as the leading character.

#### Revision History
Earlier versions of SHIELDS-PTM used binary inputs, with each
component of both E and B having separate files. The grid layout
was stored in yet another file. To make it easier to ensure
consistency between files and to make the files more user-friendly
the format was changed to a single ASCII file. This is what is
described in this document.

### Input Files: Run Configuration

Three different files are required to fully specify the test particle simulation:
```ptm_parameters_rrrr.txt```, ```dist_density_rrrr.txt```, and ```dist_velocity_rrrr.txt```.
The first of these set the options for the simulation, including particle mass and charge,
direction of time integration, integrator, etc. The latter two specify the distribution of
test particles.

#### Simulation Parameters

#### Test Particle Distribution


### Output Files: Trajectories

### Output Files: Flux Map
