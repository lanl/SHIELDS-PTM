import bisect
import os
import glob
import re
import argparse
import datetime as dt
import numpy as np
import spacepy.datamodel as dm
import spacepy.time as spt
import spacepy.coordinates as spc

import gpstools as gpt

gpspath = '/n/space_data/LANL_GPS/Particle_Data/processed_ascii'


def findDataFile(sat, targ_time, version='1.08'):
    """Find the data file that has the requested data"""
    def wkstr_to_date(wkstr):
        # Convert datestring from GPS file title to year, month, day
        year = 2000+int(wkstr[:2])
        mon = int(wkstr[2:4])
        day = int(wkstr[4:6])
        return dt.datetime(year, mon, day)
    # Get all filenames for this satellite and convert week start to date
    allfiles = sorted(glob.glob(os.path.join(gpspath, sat, '*_v{}.ascii'.format(version))))
    weekstrs = [re.search('(?<=_)\d{6}', fn).group(0) for fn in allfiles]
    weekdates = [wkstr_to_date(wk) for wk in weekstrs]
    # Find where the requested date would be inserted
    fidx = bisect.bisect(weekdates, targ_time)
    # Return filename to the left of that in time
    return allfiles[fidx-1]


def getPosition(sat, targ_time, verbose=False):
    """Get GSM position of given GPS satellite at requested time
    """
    gpsfile = findDataFile(sat, targ_time)
    data = dm.readJSONheadedASCII(gpsfile)
    # Get time (this is GPS time, not UTC, even though we're treating it that way)
    utc = spt.doy2date(data['year'].astype(int), data['decimal_day'], dtobj=True, flAns=True)
    # Find time nearest to requested
    idx = bisect.bisect(utc, targ_time)
    d1 = np.abs((utc[idx-1]-targ_time).total_seconds())
    d2 = np.abs((utc[idx]-targ_time).total_seconds())
    idx = idx-1 if d1 <= d2 else idx
    # Now get the position and convert to GSM
    lon = data['Geographic_Longitude'][idx]
    lat = data['Geographic_Latitude'][idx]
    rad = data['Rad_Re'][idx]
    cc = spc.Coords([rad, lat, lon], 'GEO', 'sph')
    cc.ticks = spt.Ticktock(utc[idx])
    pos = cc.convert('GSM', 'car')

    # Return position as a strnig that can be used in makeRun script
    posargs = '-p {0:.3f} {1:.3f} {2:.3f}'.format(pos.x[0], pos.y[0], pos.z[0])
    if verbose:
        print(sat, data['L_shell'][idx])
        print(posargs)

    return posargs


def getSpectrum(sat, targ_time):
    """Get energy versus flux of given GPS satellite at requested time
    """
    gpsfile = findDataFile(sat, targ_time)
    data = dm.readJSONheadedASCII(gpsfile)
    # Get UTC time
    data['UTC'] = gpt.computeTime(data['year'], data['decimal_day']).UTC
    # Find time nearest to requested
    idx = bisect.bisect(data['UTC'], targ_time)
    d1 = np.abs((data['UTC'][idx-1]-targ_time).total_seconds())
    d2 = np.abs((data['UTC'][idx]-targ_time).total_seconds())
    idx = idx-1 if d1 <= d2 else idx
    # Get energy spectrum
    energy = data['proton_flux_fit_energy'][idx]
    flux = data['proton_flux_fit'][idx]
    return {'energy': energy, 'flux': flux}


def get_response(svn):
    response = gpt.response.read_response(svn, species='proton')
    return response



if __name__ == '__main__':
    # Set up a basic argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sat', dest='sat', type=int,
                        help='Satellite vehicle number (as integer)')
    parser.add_argument('-t', '--time', dest='time', type=int, nargs='+',
                        help='Date time as whitespace-delimited list,'
                        + ' e.g., "-t 2017 9 3 13 15"')
    parser.add_argument('-l', dest='add_l', action='store_true')
    # To use the script directly with makeRun.py
    # Ensure that Spacepy is not returning the support notice
    # the call makeRun like:
    # python makeRun.py -o ~/test_PTM $(python scripts/gps_position.py -s 54 -t 2017 9 6 13 15)

    opt = parser.parse_args()
    if opt.time is None:
        targ_time = dt.datetime(2017, 9, 6, 13, 15)
    else:
        targ_time = dt.datetime(*opt.time)

    if opt.sat is None:
        sats = ['ns{}'.format(sn) for sn in range(54, 74)]
        targ_time = dt.datetime(2017,9,6,13,15)
        for sat in sats:
            getPosition(sat, targ_time, verbose=True)
    else:
        posargs = getPosition('ns{}'.format(opt.sat), targ_time, verbose=opt.add_l)
        if not opt.add_l: print(posargs)
