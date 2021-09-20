# Converts a 3.7G PTM data/output into a ~250M HDF5 file
# Takes roughly 4 minutes to process a single run
import sys
import glob
import spacepy.datamodel as dm
from ptm_python import ptm_preprocessing
from ptm_python import ptm_postprocessing
from ptm_python import ptm_tools

# Hardcoded path shortcuts
ptm_input  = "/ptm_input/"
ptm_output = "/ptm_output/"
ptm_data   = "/ptm_data/"

# Input
d = sys.argv[1]
run_dir = d.split("/")[-1]

# Files
map_files   = glob.glob(d + ptm_output + "map_*")
trace_files = glob.glob(d + ptm_output + "ptm_*")
field_files = glob.glob(d + ptm_data   + "ptm_fields_*")

data = dm.SpaceData()
# Parse rundir for this information?
# Or is it stored somewhere accessible?
data.attrs["time"]  = ""
data.attrs["model"] = ""
data.attrs["sat"]   = ""

# Process the field data
i = 0 # How to get i?
field = dm.SpaceData()
field_file = field_files[i]
field_data = ptm_preprocessing.PTMfields.from_file(field_file) # comes out as python class, which is kind of annoying
field_keys = ["x", "y", "z", "bx", "by", "bz", "ex", "ey", "ez", "nx", "ny", "nz"]
for key in field_keys:
    field[key] = dm.dmarray(field_data.__getattribute__(key))

# Process the map data
maps = dm.SpaceData()
map_data = ptm_tools.parse_map_file(map_files) # comes out as a regular dictionary
for k, n in map_data.items():
    maps[k] = dm.dmarray(n)

# Process the trace data
trace = dm.SpaceData()
for fname in trace_files:
    trace_data = ptm_tools.parse_trajectory_file(fname) # comes out as a regular dictionary
    for k, n in trace_data.items(): # dictionary keys are particle numbers, can we concatenate dictionaries together?
        trace[k] = dm.dmarray(n).T # why do I want them transposed?

data["Field"] = field
data["Map"]   = maps
data["Trace"] = trace

data.toHDF5(d + "/PTM.h5")
