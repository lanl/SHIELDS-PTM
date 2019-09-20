-------------------------------------------------------
| ANALYSIS OF DATA FROM THE SHIELDS-PTM PARTICLE CODE |
-------------------------------------------------------
|                        |
| Jesse Woodroffe, ISR-1 |
| jwoodroffe@lanl.gov    |
|                        |
--------------------------


I. Background

The SHIELDS-PTM particle tracing code solves the equations of charged particle
motion in either the guiding center drift or full-particle (Lorentz) orbit 
representation.

The SHIELDS-PTM simulation is configured using three text files:

*ptm_parameters_xxxx.txt
*dist_density_xxxx.txt
*dist_velocity_xxxx.txt

where xxxx stands for a four-digit, left-zero-padded integer (such as 0001)

The SHIELDS-PTM simulation outputs (up to) two different types of output files:

*ptm_xxxx.dat
*map_xxxx.dat

where xxxx stands for the same four digit zero-padded integer as with the input file

The "ptm" output file contains trajectories (at a cadence spencified in ptm_parameters_xxxx.txt)
for each particle in the simulation. For each time step that is output into this file, there are
seven quantities that are recorded:

TIME	X-Posiition	Y-Position	Z-Position	Vperp	Vpara	Energy  PitchAngle

When SHIELDS-PTM is run in parallel, particles are not output in any particular order. Instead,
a particle's full time history is output as soon as it is available (i.e. as soon as the particle
has finished integration). See the description in file_format.txt for a complete description of
how these files are arranged.

Because of the inherent disorder of the ptm_xxxx.dat file, it is not convenient to use this
file for Liouville tracing estimates of the distribution, which is one of the primary reasons
that SHIELDS-PTM was developed in the first place. However, the map_xxxx.dat file is designed
for just this purpose, and so 

II. Units



III. Derived Quantities




IV. 
