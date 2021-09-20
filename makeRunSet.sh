export OMP_NUM_THREADS=20

# Directory to save the runs to
runsdir='/projects/lanl/Carrington/ModelOutputs/PTM/test'

# Directory to save the fields to
fieldsdir='/projects/lanl/Carrington/ModelInputs/PTMFields/test'

# Directory for PTM
ptmdir=`pwd`

# Date to run
date='20170907_050000'

# Size of time steps (seconds)
deltat=300

# Number of time steps
nt=1

# Satellites to simulate
sats='57' # 62 68 71 72'

# Models to simulate
models='CDIP' # T96 T01S T02 TS04'

# Step 1: Get variables to make stuff easier
ptmtools=$ptmdir/tools
ptmpython=$ptmdir/ptm_python
ptmscripts=$ptmdir/scripts
day=${date:0:8}
tim=${date:9:8}
xt=$(expr $deltat \* 100 \/ 60)
dates=$date
for ((i=1; i<nt; i++)); do
  y=$(expr $xt \* $i)   # This method of getting times only works if you don't
  d=$(expr $tim \+ $y)  # go over the hour mark. And if the deltat is divisible
  printf -v t "%06d" $d # by 60.
  dates="${dates} ${day}_$t"
done

# Step 2: Create the fields for the different models
mkdir -p $fieldsdir
for model in $models; do
  dir=$fieldsdir/$date-$model
  mkdir -p $dir
  rm $dir/*
  echo $ptmtools/ptm_lanlgeomag -d $day -t $tim -e $model -n $nt -l $deltat $dir
  $ptmtools/ptm_lanlgeomag -d $day -t $tim -e $model -n $nt -l $deltat $dir
done

# Step 3: Make PTM run directories
mkdir -p $runsdir
for d in $dates; do
  for model in $models; do
    for sat in $sats; do
      dir=$runsdir/$d-$model-ns$sat
      mkdir -p $dir
      mkdir -p $dir/ptm_output
      mkdir -p $dir/ptm_input
      mkdir -p $dir/ptm_data
      rm $dir/ptm
      rm $dir/ptm_data/*
      ln -s $ptmdir/ptm $dir/.
      ln -s $fieldsdir/$date-$model/* $dir/ptm_data/.
      echo python3 $ptmscripts/gps_position.py --time ${d:0:4} ${d:4:2} ${d:6:2} ${d:9:2} ${d:11:2} --sat $sat -l
      python3 $ptmscripts/gps_position.py --time ${d:0:4} ${d:4:2} ${d:6:2} ${d:9:2} ${d:11:2} --sat $sat -l > $dir/gps_position.log
      tag=$(tail -n 1 $dir/gps_position.log)
      seconds=$(expr $i \* $deltat)
      echo python3 $ptmdir/makeRun.py -r 1 -o $dir -s $seconds $tag
      python3 $ptmdir/makeRun.py -r 1 -o $dir -s $seconds $tag
    done
  done
done

# Step 4: Run PTM
for d in $dates; do
  for model in $models; do
    for sat in $sats; do
      dir=$runsdir/$d-$model-ns$sat
      cd $dir
      bash job_script.sh
      cd $ptmdir
      echo python3 $ptmscripts/PTMtoHDF5.py $dir
      python3 $ptmscripts/PTMtoHDF5.py $dir
    done
  done
done
