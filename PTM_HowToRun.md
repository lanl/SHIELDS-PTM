# Quickstart guide for PTM

Any shell commands in-text are given in backticks (or rendered in a code environment).
Further details of how PTM (Particle Tracing Model) works are given in `PTM_doc`

## A. Prepare PTM input files

   Python3 scripts to prepare PTM input files should be in directory scripts/.
   Create there directories `ptm_data/` and `ptm_input/` and keep there the .mhd magnetic-field--model files

   1. `module load python/3.6-anaconda-5.0.1` or any of the python3/anaconda modules 
   2. `python ptm_tec_interp.py` to create PTM input data files in `ptm_data/` by interpolating the .mhd files 
   3.  `python` to create PTM input parameter files in `ptm_input/`

      ```
      import ptm_input
      ptm_input.ptm_input_creator().create_input_files()
      ```

      or 

      ```
      ptm_input.ptm_input_creator().get_interactive_input()
      ```

      if you want to select all input parameters interactively

## B. Running PTM

### Single process

   1. load module gfortran, if needed
   2. create a directory `ptm_input/` containing all 3 parameter files created before by `PTM_INPUT.py`
   3. create a directory `ptm_data/` containing all 4+3*2*n data files created before by `PTM_TEC_INTERP.py`
   Note: PTM needs B&E input files at least two epochs (n >=2). If you have only one input .mhd file,
         then duplicate the only interpolated (b/e)(x,y,z)3d__0001.dat (n=1) files and call them 
         (b/e)(x,y,z)_0002.dat (n=2). In this magnetostatic set-up, PTM will do only spatial interpolations
         to calculate B&E at current particle location.
   4. `make all` to create executable ptm
   5. `./ptm N` select N (=runid). 
   Note: PTM will read input files in directory `ptm_input/` ending in 000N.txt
                  write output file `ptm_output/ptm_000N.dat`

### Multiprocessing
   1. Make a run directory on scratch, e.g.,

        setenv PTM_RUNDIR /net/scratch4/${USER}/PTMrun
        mkdir /net/scratch4/${USER}/PTMrun

   2. Copy PTM input data, parameter files, executable, to run directory

        cp -r ptm_input ${PTM_RUNDIR}
        cp -r ptm_data ${PTM_RUNDIR}
        cp ptm ${PTM_RUNDIR}

   3. Write SLURM job script (say, `job_script.sh`) to launch jobs in parallel using GNU parallel. E.g.

        #!/bin/tcsh
        #SBATCH --time=1:00:00
        #SBATCH --nodes=2
        #SBATCH --qos=standard
        #SBATCH --account=w19_gic
        #SBATCH --job-name=PTM_test
        #SBATCH --no-requeue
        #SBATCH -o slurm_%x_%j.out
        #SBATCH -e slurm_%x_%j.err
        #SBATCH --mail-user=username@lanl.gov
        #SBATCH --mail-type=BEGIN,END,FAIL

        module purge
        module load gcc/6.4.0

        parallel --delay 0.2 -j $SLURM_NTASKS \
                 "srun -N 1 -n 1 ./ptm" ::: \
                 1 2 3 4 5 #This final line is a whitespace-separated list of input arguments for PTM

   5. Save the SLURM job submission script in `${PTM_RUNDIR}`.

   6. Submit the job script

        sbatch job_script.sh

## C. Post-process PTM output

   PTM output is post-processed using Python 3 scripts.

     python
     >>> import ptm_postprocessing
     >>> p=ptm_postprocessing.ptm_postprocessor()
     >>> p.process_run(N) (N must match that of existing PTM output file map_000N.dat)

   to process Flux-Map output file `map_000N.dat` created with PTM for idist=3,4 in `dist_velocity_000N.txt`
