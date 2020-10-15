import glob
import sys
import subprocess

cadence = 300
newtime = sys.argv[1]
lostring = r's/^.*\(tlo.*$\)/{}      \1/'.format(int(newtime)-300)
histring = r's/^.*\(thi.*$\)/{}      \1/'.format(newtime)
#dtstring = r's/^.*\(dtin.*$\)/{}      \1/'.format(cadence)
#ilstring = r's/^.*\(ilast.*$\)/{}      \1/'.format(2)
#ntstring = r's/^.*\(ntot.*$\)/{}      \1/'.format(2)

flist = glob.glob('ptm_parameters_*.txt')
for fname in flist:
    subprocess.call(['sed', '-i', lostring, fname])
    subprocess.call(['sed', '-i', histring, fname])
    #subprocess.call(['sed', '-i', dtstring, fname])
    #subprocess.call(['sed', '-i', ilstring, fname])
    #subprocess.call(['sed', '-i', ntstring, fname])
