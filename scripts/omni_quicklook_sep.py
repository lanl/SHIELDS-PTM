import sys
import os
import glob
import itertools as it
import numpy as np
import spacepy.toolbox as tb
import matplotlib.pyplot as plt

# Because ptm_python is currently not an installed/installable,
# module,  scripts should set the location so it can be
# imported
sys.path.insert(1, os.path.dirname(os.path.abspath('../ptm_python')))
from ptm_python import ptm_tools as ptt
from ptm_python import ptm_postprocessing as post

def calculate_omni(fns):
    fluxmap = ptt.parse_map_file(fns)
    pp = post.ptm_postprocessor()
    
    pp.set_source_parameters(n_dens=5e-6, e_char=500.0, kappa=4.0, mass=1847.0)
    # n_dens      float       optional, number density at source region in cm-3
    # e_char      float       optional, characteristic energy of distribution in keV
    # kappa       float       optional, spectral index of kappa distribution
    # mass        float       optional, particle mass in multiples of electron mass
    
    diff_flux = pp.calculate_flux(fluxmap)
    nener, nalph = diff_flux.shape
    for idx, jdx in it.product(range(nener), range(nalph)):
        if (np.abs(fluxmap['final_x'][idx, jdx]) >= 15).any():
            continue
        else:
            diff_flux[idx, jdx] = 0.
    
    omni = pp.calculate_omnidirectional_flux(fluxmap['angles'], diff_flux)

    return omni, fluxmap

def plot_omni(omni, fluxmap):
    fig = plt.figure()
    ax0 = fig.add_axes([0.15, 0.2, 0.78, 0.6])
    en_mev = fluxmap['energies']/1e3
    ax0.loglog(en_mev, omni/1e3)
    enlo, enhi = en_mev[0], en_mev[-1]
    ax0.set_xlabel('Energy [MeV]')
    ax0.set_ylabel('Diff. Flux [per MeV]')
    ax0.set_title(fluxmap.attrs['position'])
    
    # find energies where access is forbidden
    if False:  # binary - any access at E or no access
        allow = np.array([True if (np.abs(fluxmap['final_x'][n])==15).any()
                          else False for n in range(len(fluxmap['energies']))])
    else:  # fractional
        allow = np.array([np.sum(np.abs(fluxmap['final_x'][n])==15)
                          for n in range(len(fluxmap['energies']))], dtype=np.float)
        allow /= len(fluxmap['angles'])
    def cutoffs(allow):
        cdict = dict()
        not_forbidden = np.nonzero(allow)[0]
        idx_low = not_forbidden[0]
        ec_low = en_mev[idx_low]
        idx_high = not_forbidden[-1]
        ec_high = en_mev[idx_high]
        prot_low = ptt.Proton(ec_low)
        prot_high = ptt.Proton(ec_high)
    
        Nforbid = np.sum(allow*len(fluxmap['angles']))
        Ntot = len(allow)*len(fluxmap['angles'])
        r_low = prot_low.getRigidity()
        r_high = prot_high.getRigidity()
        r_effect = r_low + (r_high-r_low) * Nforbid/Ntot
        cdict['ec_low'] = ec_low
        cdict['r_low'] = r_low
        cdict['ec_high'] = ec_high
        cdict['r_high'] = r_high
        cdict['r_eff'] = r_effect
        cdict['ec_eff'] = ptt.Proton.fromRigidity(r_effect).energy
        return cdict
    
    cdict = cutoffs(allow)
    #ax0.axvline(x=cdict['ec_low'], linestyle=':', color='k',
    #            label='E$_{c}^{low}$ = ' + '{0:.2f} MeV'.format(cdict['ec_low']))
    #ax0.axvline(x=cdict['ec_high'], linestyle=':', color='k',
    #            label='E$_{c}^{high}$ = ' + '{0:.2f} MeV'.format(cdict['ec_high']))
    ax0.axvline(x=cdict['ec_eff'], linestyle='--', color='b',
                label='E$_{c}^{eff}$ = ' + '{0:.2f} MeV'.format(cdict['ec_eff'])
                + '\nR$_C$ = ' + '{0:.2f} GV'.format(cdict['r_eff']))
    ax0.legend()
    # a vertical barcode
    ax1 = fig.add_axes([0.15, 0.05, 0.78, 0.033])
    ax1.set_axis_off()
    #barprops = dict(aspect='auto', cmap='binary', interpolation='nearest')
    barprops = dict(cmap='gray')
    binedges = tb.bin_center_to_edges(en_mev)
    ax1.pcolormesh(np.log10(binedges), np.array([0, 1]),
                   allow[np.newaxis, ...], **barprops)
    enlo, enhi = binedges[0], binedges[-1]
    ax0.set_xlim([enlo, enhi])
    ax1.set_xlim(np.log10([enlo, enhi]))
    
    plt.show()

if __name__ == '__main__':
    fns = glob.glob('/net/toaster/u/smorley/test_PTM/ptm_output/map_00*.dat')
    omni, fluxmap = calculate_omni(fns)
    plot_omni(omni, fluxmap)
