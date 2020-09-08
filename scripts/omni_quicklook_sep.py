import sys
import os
import glob
import copy
import argparse
import datetime as dt
import itertools as it
import numpy as np
from scipy import interpolate
import spacepy.toolbox as tb
import matplotlib.pyplot as plt

from ptm_python import ptm_tools as ptt
from ptm_python import ptm_postprocessing as post

import gps_position


class CXD(post.ptm_postprocessor):
    """CXD Instrument subclassed from PTM post-processor
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cutoff_info = dict()
        self.instrument_info = dict()


    def set_source(self, source_type='kappa', params={}, empflux={}):
        '''Set source flux spectrum [cm^-1 s^-1 keV^-1 sr^-1]

        If empirical, pass in a dict of energy and flux values
        empflux = {'flux': [values], 'energy': [values]}
        '''
        if source_type == 'empirical':
            self.get_flux = interpolate.interp1d(empflux['energy'], empflux['flux'],
                                                 kind='linear', axis=-1, copy=True,
                                                 bounds_error=None, fill_value='extrapolate',
                                                 assume_sorted=False)
            self.source_type = source_type
        else:
            super().set_source(source_type=source_type, params=params)


    def calculate_omni(self, fluxmap, fov=False, initialE=False, source=None):
        '''Source flux spectrum is per steradian. Output is per steradian.
        '''
        if source is not None:
            self.set_source(source_type='empirical', empflux=source)
        else:
            self.set_source(source_type='kaprel', params=dict(density=5e-6, energy=752.0, kappa=5.0, mass=1847.0))
            # n_dens      float       optional, number density at source region in cm-3
            # e_char      float       optional, characteristic energy of distribution in keV
            # kappa       float       optional, spectral index of kappa distribution
            # mass        float       optional, particle mass in multiples of electron mass
        if initialE:
            fluxmap['final_E'] = fluxmap['init_E']
        diff_flux = self.map_flux(fluxmap)
        if not initialE:
            nener, nalph = diff_flux.shape
            for idx, jdx in it.product(range(nener), range(nalph)):
                if np.linalg.norm(fluxmap['final_x'][idx, jdx]) >= 14.99:
                    continue
                else:
                    diff_flux[idx, jdx] = 0.
        per_ster = 4*np.pi
        if fov:
            self.instrFOV(fluxmap, fov=fov)
            diff_flux[self.instrument_info['angle_mask']] = 0.0
            per_ster = self.instrument_info['norm_per_steradian']
        omni = self.get_omni_flux(fluxmap['angles'], diff_flux)
        # omni, at this point, is integrated over look direction
        # so we divide by a normalization factor to get our flux per steradian
        omni /= per_ster

        return omni


    def instrFOV(self, fluxmap, fov=90, dir='nadir', verbose=False):
        # get the instrument field-of-view axis
        # CXD is nadir pointing, so FOW axis is inverse of position vector
        posvec = fluxmap.attrs['position']
        posvec = posvec/np.linalg.norm(posvec)
        fovaxis = -1*posvec
        fovaxis = fovaxis/np.linalg.norm(fovaxis)
        nadir = fovaxis.copy()
        if dir.lower() == 'east':
            fovaxis = np.cross(fovaxis, [0, 0, 1])
            fovaxis = fovaxis/np.linalg.norm(fovaxis)
        elif dir.lower() == 'west':
            fovaxis *= -1
        elif dir.lower() == 'nadir':
            pass
        else:
            raise NotImplementedError('Requested direction ({}) not defined.'.format(dir))
        if verbose:
            print("Position: {:0.3f},{:0.3f},{:0.3f}".format(*posvec.tolist()))
            print("Nadir: {:0.3f},{:0.3f},{:0.3f}".format(*nadir.tolist()))
            print("FOV: {:0.3f},{:0.3f},{:0.3f}".format(*fovaxis.tolist()))
        # get angle between inverse of FOV axis and initial particle velocity
        # (because velocity needs to point in to the detector aperture)
        fiv = fluxmap['init_v']
        angles = np.zeros((fiv.shape[0], fiv.shape[1]))
        for idx, jdx in it.product(range(fiv.shape[0]), range(fiv.shape[1])):
            angles[idx, jdx] = np.rad2deg(np.arccos(np.dot(-1*fovaxis, fluxmap['init_v'][idx, jdx]
                                                           /np.linalg.norm(fluxmap['init_v'][idx, jdx]))))
        self.instrument_info['angles'] = angles
        self.instrument_info['angle_mask'] = angles < fov
        self.instrument_info['FOV_degrees'] = fov
        # self.instrument_info['norm_per_steradian'] = 2*np.pi*(1-np.cos(np.deg2rad(fov)))
        # Normalization is borked because omnidirectional flux integration assumes gyrotropy
        if fov == 90:
            self.instrument_info['norm_per_steradian'] = 2*np.pi
        else:
            NotImplementedError


    def cutoffs(self, fluxmap, addTo=None, verbose=True, **kwargs):
        # find energies where access is forbidden
        if False:  # binary - any access at E or no access
            allow = np.array([True if (np.abs(fluxmap['final_x'][n])==15).any()
                              else False for n in range(len(fluxmap['energies']))])
        else:  # fractional
            allow = np.array([np.sum(np.linalg.norm(fluxmap['final_x'][n], axis=-1) >= 14.99)
                              for n in range(len(fluxmap['energies']))], dtype=np.float)
            allow /= len(fluxmap['angles'])
    
        en_mev = fluxmap['energies']/1e3
        cdict = dict()
        # Earth (1.1Re to include atmosphere) subtends ~0.225 ster
        # (solid angle = 2*pi*(1-cos(theta)), where theta is angle of inverse
        # position vector to limb of Earth)
        # 0.225 ster/4*pi = 0.179
        not_forbidden = np.nonzero(allow)[0]
        full_access = allow >= 0.821  # 1-0.179 = 0.821
        idx_low = not_forbidden[0]
        ec_low = en_mev[idx_low]
        # upper cutoff is where all values are "full" transmission
        # Here "full" accounts for solid angle of Earth
        try:
            idx_high = find_runs(full_access)[-1][0]
        except IndexError:
            idx_high = -1
        ec_high = en_mev[idx_high]
        if verbose:
            print('Ec_low = {}'.format(ec_low))
            print('Ec_high = {}'.format(ec_high))
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
        cdict['allow'] = allow
        if addTo is not None:
            if 'linestyle' not in kwargs:
                kwargs['linestyle'] = '--'
            if 'color' not in kwargs:
                kwargs['color'] = 'dimgrey'
            labeltext = 'E$_{c}^{eff}$ = ' + '{0:.2f} MeV'.format(cdict['ec_eff'])\
                        + '\nR$_C$ = ' + '{0:.2f} GV'.format(cdict['r_eff'])
            if 'label_pre' in kwargs:
                labeltext = kwargs['label_pre'] + '\n' + labeltext
            addTo.axvline(x=cdict['ec_eff'], linestyle=kwargs['linestyle'], color=kwargs['color'],
                          label=labeltext)
        self.cutoff_info = cdict
        return cdict, allow


def plot_omni(instr, omni, fluxmap, label='Omni'):
    fig = plt.figure()
    ax0 = fig.add_axes([0.15, 0.2, 0.78, 0.6])
    en_mev = fluxmap['energies']/1e3
    ax0.loglog(en_mev, omni*1e3, label=label)
    enlo, enhi = en_mev[0], en_mev[-1]
    ax0.set_xlabel('Energy [MeV]')
    ax0.set_ylabel('Diff. Flux [per MeV]')
    ax0.set_title(fluxmap.attrs['position'])

    cdict, allow = instr.cutoffs(fluxmap, addTo=ax0, linestyle='--', color='b')
    #ax0.axvline(x=cdict['ec_low'], linestyle=':', color='k',
    #            label='E$_{c}^{low}$ = ' + '{0:.2f} MeV'.format(cdict['ec_low']))
    #ax0.axvline(x=cdict['ec_high'], linestyle=':', color='k',
    #            label='E$_{c}^{high}$ = ' + '{0:.2f} MeV'.format(cdict['ec_high']))

    # a horizontal barcode
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
    return fig, [ax0, ax1]



def find_runs(invec, value=True):
    """Find starts of contiguous blocks"""
    isvalue = np.concatenate(([0], np.equal(invec, value).view(np.int8), [0]))
    absdiff = np.abs(np.diff(isvalue))
    # runs start and end where absdiff is 1
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges


def add_extra_omni(ax, omni, fluxmap, **kwargs):
    ax.loglog(fluxmap['energies']/1e3, omni*1e3, **kwargs)



if __name__ == '__main__':
    import bisect
    import gpstools as gpt

    # Set up a basic argument parser
    parser = argparse.ArgumentParser()
    # Add a positional argument for the ptm_output folder
    parser.add_argument('dir')
    opt = parser.parse_args()
    # Grab the simulation output and set up post-processor
    fns = glob.glob(os.path.join(opt.dir, 'map_*.dat'))
    cxd = CXD()
    fluxmap = ptt.parse_map_file(fns)
    fluxmap_inst = copy.deepcopy(fluxmap)
    fluxmap_source = copy.deepcopy(fluxmap)
    # Get source spectrum from minimally-shielded CXD
    targ = dt.datetime(2017, 9, 6, 13, 5)
    obsdict = gps_position.getSpectrum('ns61', targ)
    obsdict['energy'] = obsdict['energy']*1000  # Convert energy from MeV to keV
    obsdict['flux'] = obsdict['flux']/1000  # Convert flux from per MeV to per keV
    fluxlim = 60  # MeV - lower limit for integration
    # Get target spectrum
    obs = gpt.loadCXDascii(72, '17090*')
    idx = bisect.bisect_left(obs['UTC'], targ)
    obs_gt = interpolate.BSpline(*(interpolate.splrep(obs['proton_flux_fit_energy'][idx],
                                   obs['proton_flux_fit'][idx]))).integrate(np.min(fluxlim),np.max(800))
    # Calculate fluxes and plot
    omni = cxd.calculate_omni(fluxmap, source=obsdict)
    omni_inst = cxd.calculate_omni(fluxmap_inst, source=obsdict, fov=90)
    omn_gt = interpolate.BSpline(*(interpolate.splrep(fluxmap['energies'], omni_inst)
                                 )).integrate(np.min(fluxlim*1000),np.max(800000))
    fig, axes = plot_omni(cxd, omni_inst, fluxmap, label='ns72 predicted')
    omni_source = cxd.calculate_omni(fluxmap_source, source=obsdict, initialE=True)
    print(">{}MeV: obs={}; pred={}".format(fluxlim, obs_gt, omn_gt))
    # add_extra_omni(axes[0], omni_inst, fluxmap_inst, label='CXD', color='orange')
    add_extra_omni(axes[0], omni_source, fluxmap_source, label='Source (@ ns61)', color='green')
    ## Plot initial spectrum
    # axes[0].plot(obsdict['energy']/1000, obsdict['flux']*1000, 'k-.')
    ## Plot expected, everything is omni flux _per steradian_
    # axes[0].plot(obs['proton_flux_fit_energy'][idx], obs['proton_flux_fit'][idx], 'y:', label='ns72')
    axes[0].legend()
    ylims = axes[0].get_ylim()
    cdict = cxd.cutoff_info
    axes[0].plot(cdict['ec_low'], ylims[0], marker='^', mec='k', mfc='silver', clip_on=False)
    axes[0].plot(cdict['ec_high'], ylims[0], marker='^', mec='k', mfc='grey', clip_on=False)
    axes[0].set_ylim(ylims)
    plt.show()
