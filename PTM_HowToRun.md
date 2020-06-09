# Quickstart guide for PTM

Any shell commands in-text are given in backticks (or rendered in a code environment).
Further details of how PTM (Particle Tracing Model) works are given in `PTM_DOC`

## A. Prepare PTM input files

   Magnetic and electric field input files are required by PTM for tracing particles.
   For the purposes of this quickstart we assume that these files will be generated from
   Space Weather Modeling Framework output.

   Preparing the PTM input files is done using the `ptm_preprocessing.py` Python script.
   Supported SWMF filetypes include the various `idl` file formats, as well as the ASCII
   Tecplot format. The `idl` binary files are preferred due to their smaller size.
   Either `MHD` or `FUL` outputs should be used as these provide the required information
   for PTM.

   To use the `ptm_python` tools, first ensure that you have a Python 3 environment with
   numpy, scipy, matplotlib and spacepy installed. On institutional computing, this is
   accomplished with modules, e.g.:

   1. `module load python/3.6-anaconda-5.0.1`
   2. `python ptm_preprocessing.py` to create PTM input data files in `ptm_data/` by interpolating the SWMF 3d output files 
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

      A testcase for SEP access is provided. Simply use the `makeRun.py` script and
      this will generate the required input files. Help can be obtained by using the
      `--help` option, e.g., `python makeRun.py --help`

## B. Running PTM

### Single process

   1. Ensure your environment has the same dependencies as at build time. E.g., `module load gcc/8.1.0`
   2. create a directory `ptm_input/` containing all 3 parameter files created before by `ptm_input.py`
   3. create a directory `ptm_data/` containing all 4+3*2*n data files created before by `ptm_preprocessing.py`
   Note: PTM needs B&E input files at least two epochs (n >=2). If you have only one input SWMF file,
         then ptm_preprocessing will process it twice, making two sets of output files. This allows
         in-time interpolation for a static field.
   4. `make all` to create executable ptm
   5. `./ptm N` select N (=runid). 
   Note: PTM will read input files in directory `ptm_input/` ending in ????.txt
                  write output file `ptm_output/ptm_????.dat`

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
     >>> p = ptm_postprocessing.ptm_postprocessor()
     >>> p.process_run(N)  # (N must match that of existing PTM output file map_000N.dat)

   to process Flux-Map output file `map_000N.dat` created with PTM for idist=3,4 in `dist_velocity_000N.txt`


# Sample Workflow
   First we preprocess 45 minutes of SWMF output:
    (default) [ptm_python]$ python ptm_preprocessing.py -t 900 -o ~/projects/SEP_ER/ns73_1 /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-13*.out
    start reading /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-130000-001.out
    done reading file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-130000-001.out
    done interpolating file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-130000-001.out
    start reading /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-131500-020.out
    done reading file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-131500-020.out
    done interpolating file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-131500-020.out
    start reading /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-133000-001.out
    done reading file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-133000-001.out
    done interpolating file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-133000-001.out
    start reading /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-134500-009.out
    done reading file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-134500-009.out
    done interpolating file /n/projects/lanl/Carrington/ModelOutputs/SWMF_Results/Sept2017_hi_SEP1/GM/3d__mhd_1_e20170906-134500-009.out

   Then we set up runs 10 minutes apart using positions along a GPS orbit (2 shown)
    (default) [SHIELDS-PTM]$ mkdir ~/projects/SEP_ER/ns72/20170906133500
    (default) [SHIELDS-PTM]$ ln -s ~/projects/SEP_ER/Sept2017_SEP1/ ~/projects/SEP_ER/ns72/20170906133500/ptm_data
    (default) [SHIELDS-PTM]$ python scripts/gps_position.py --time 2017 9 6 13 35 --sat 72 -l
    ns72 5.022402
    -p 3.356 2.246 1.076
    (default) [SHIELDS-PTM]$ python makeRun.py -r 1 -o ~/projects/SEP_ER/ns72/20170906133500 -s 2100 -p 3.356 2.246 1.076
    Run jobs using `source job_script.sh` from the run directory

   ...

    (default) [SHIELDS-PTM]$ mkdir ~/projects/SEP_ER/ns72/20170906130500
    (default) [SHIELDS-PTM]$ ln -s ~/projects/SEP_ER/Sept2017_SEP1/ ~/projects/SEP_ER/ns72/20170906130500/ptm_data
    (default) [SHIELDS-PTM]$ python scripts/gps_position.py --time 2017 9 6 13 05 --sat 72 -l
    ns72 6.170052
    -p 3.408 1.597 1.817
    (default) [SHIELDS-PTM]$ python makeRun.py -r 1 -o ~/projects/SEP_ER/ns72/20170906130500 -s 300 -p 3.408 1.597 1.817
    Run jobs using `source job_script.sh` from the run directory

