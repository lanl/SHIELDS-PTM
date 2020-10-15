import sys
import os
import glob
import copy
import argparse
import datetime as dt
import itertools as it
import numpy as np
from scipy import interpolate, integrate
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


    def calculate_diff(self, fluxmap, initialE=False, source=None):
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
        fluxmap['diff_flux'] = self.map_flux(fluxmap)
        if not initialE:
            nener, nalph = fluxmap['diff_flux'].shape
            self.instrument_info['allowed'] = np.full(fluxmap['diff_flux'].shape, False, dtype=np.bool_)
            for idx, jdx in it.product(range(nener), range(nalph)):
                if np.linalg.norm(fluxmap['final_x'][idx, jdx]) >= 14.99:
                    self.instrument_info['allowed'][idx, jdx] = True
                    continue
                else:
                    fluxmap['diff_flux'][idx, jdx] = 0.
        return fluxmap


    def calculate_omni(self, fluxmap, fov=False, initialE=False, source=None,
                       from_look=True, dir='nadir'):
        '''Source flux spectrum is per steradian. Output is per steradian.
        '''
        fluxmap = self.calculate_diff(fluxmap, initialE=initialE, source=source)
        per_ster = 4*np.pi
        if fov:
            self.instrFOV(fluxmap, fov=fov, dir=dir)
            fluxmap['diff_flux'][self.instrument_info['angle_mask']] = 0.0
            per_ster = self.instrument_info['norm_per_steradian']
        if from_look:
            nener = len(fluxmap['energies'])
            fluxmap['gridded_j'] = np.zeros([nener, 200, 400])
            grid_theta, grid_phi = np.mgrid[0:np.pi:200j, -np.pi:np.pi:400j]
            fluxmap['grid_theta'] = grid_theta[:, 0]
            fluxmap['grid_phi'] = grid_phi[0, :]
            omni = np.zeros(nener)
            theta = np.arctan2(fluxmap['init_v'][:,:,0], fluxmap['init_v'][:,:,1])
            phi = np.arccos(fluxmap['init_v'][:,:,2]/np.linalg.norm(fluxmap['init_v'][:,:,:], axis=-1))  
            for eidx in range(nener):
                gridded_j = interpolate.griddata(np.vstack([phi[eidx,:], theta[eidx,:]]).T,
                                                 fluxmap['diff_flux'][eidx,:],
                                                 (grid_theta, grid_phi),
                                                 method='nearest')
                fluxmap['gridded_j'][eidx] = gridded_j
                # TODO: zero out flux outside the FOV... the NN interpolation lets flux bleed
                # from the valid points to outside of the FOV
                if fov:
                   center_theta = self.instrument_info['FOV_center']
                # Now do the integral
                nthet = grid_theta.shape[0]
                tmp = np.zeros(nthet)
                for ii in range(nthet):
                    tmp[ii] = integrate.simps(gridded_j[ii, :], grid_phi[0,:])
                fixomni = integrate.simps(tmp, grid_theta[:,0])
                omni[eidx] = fixomni
        else:
            omni = self.get_omni_flux(fluxmap['angles'], fluxmap['diff_flux'])
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
        self.instrument_info['FOV_center'] = fovaxis
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


    def get_response(self, svn):
        from gpstools import response
        self.instrument_info['response'] = response.read_response(svn, species='proton')


    def get_counts(self, energies, omniflux, svn=None, low=10e3, high=500e3):
        """ Get counts by integrating j(E)G(E) over E

        Energies in keV
        Flux in per keV
        (to match rest of PTM stuff)
        """
        if 'response' not in self.instrument_info:
            if svn is None:
                raise ValueError('Response function not set, please give SVN in call')
            else:
                self.get_response(svn)

        # low and high set bounds of integration
        get_j = interpolate.interp1d(energies/1e3, omniflux*1e3, fill_value=0, bounds_error=False)

        counts = dict()
        for chn in range(12,17):
            get_G = interpolate.interp1d(self.instrument_info['response']['energy'],
                                            self.instrument_info['response'][chn]/(2*np.pi),
                                            fill_value=0, bounds_error=False)

            counts[chn], err = integrate.quad(lambda e: get_j(e)*get_G(e), low/1e3, high/1e3,
                                              limit=80, epsabs=0.5)
        return counts


def plot_omni(instr, omni, fluxmap, label='Omni'):
    fig = plt.figure()
    ax0 = fig.add_axes([0.15, 0.2, 0.78, 0.6])
    en_mev = fluxmap['energies']/1e3
    ax0.loglog(en_mev, omni*1e3, label=label)
    enlo, enhi = en_mev[0], en_mev[-1]
    ax0.set_xlabel('Energy [MeV]')
    ax0.set_ylabel('Diff. Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$MeV$^{-1}$]')
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


def plot_directional_flux(energy, fluxmap1, fluxmap2):
    eidx1 = np.argmin(np.abs(fluxmap1['energies'] - energy))
    eidx2 = np.argmin(np.abs(fluxmap2['energies'] - energy))
    fig, axes = plt.subplots(nrows=1, ncols=2)
    im = axes[0].imshow(fluxmap1['gridded_j'][eidx1].T, extent=(0,np.pi,-np.pi,np.pi), origin='lower')
    clim=im.properties()['clim']
    axes[1].imshow(fluxmap2['gridded_j'][eidx2].T, extent=(0,np.pi,-np.pi,np.pi), origin='lower', clim=clim)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.5)
    fig.suptitle('{:0.2f} MeV'.format(fluxmap1['energies'][eidx1]/1e3))
    axes[0].set_ylabel('Phi')
    axes[0].set_xlabel('Theta')
    axes[1].set_xlabel('Theta')


def plot_counts_pred_obs(pred_rate, obs_rate, predlabel='Predicted', obslabel='Observed'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pred_err = 1.96*np.sqrt(pred_rate*240)/240
    obs_err = 1.96*np.sqrt(obs_rate*240)/240
    chns = range(1, len(pred_rate)+1)
    chn_labels = ['P{}'.format(cc) for cc in chns]
    ax.errorbar(chns, pred_rate, yerr=pred_err, label=predlabel,
                capsize=5, elinewidth=1.5, markeredgewidth=1.5,
                linestyle='none', marker='d', alpha=0.75)
    ax.errorbar(chns, obs_rate, yerr=obs_err, label=obslabel,
                capsize=5, elinewidth=1.5, markeredgewidth=1.5,
                linestyle='none', marker='d', alpha=0.75)
    ax.legend()
    ax.set_xticks(chns)
    ax.set_xticklabels(chn_labels)
    ax.set_xlabel('Channel')
    ax.set_ylabel('Count Rate')



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
    fluxmap_gyro = copy.deepcopy(fluxmap)
    fluxmap_inst = copy.deepcopy(fluxmap)
    fluxmap_source = copy.deepcopy(fluxmap)
    # Get source spectrum from minimally-shielded CXD
    targ = dt.datetime(2017, 9, 8, 5, 5)
    obsname = 'ns62'
    obsdict = gps_position.getSpectrum(obsname, targ)
    obsdict['energy'] = obsdict['energy']*1000  # Convert energy from MeV to keV
    obsdict['flux'] = obsdict['flux']/1000  # Convert flux from per MeV to per keV
    fluxlim = 50  # MeV - lower limit for integration
    # Get target spectrum
    satnum = 68
    obs = gpt.loadCXDascii(satnum, '17090*')
    idx = bisect.bisect_left(obs['UTC'], targ)
    obs_gt = interpolate.BSpline(*(interpolate.splrep(obs['proton_flux_fit_energy'][idx],
                                   obs['proton_flux_fit'][idx]))).integrate(np.min(fluxlim),np.max(800))
    # Calculate fluxes and plot
    omni = cxd.calculate_omni(fluxmap, source=obsdict, from_look=True)
    omni_gyro = cxd.calculate_omni(fluxmap, source=obsdict, from_look=False)
    omni_inst = cxd.calculate_omni(fluxmap_inst, source=obsdict, fov=90, dir='east', from_look=True)
    # Get counts...
    counts = cxd.get_counts(fluxmap['energies'], omni_inst, svn=satnum, low=10e3, high=500e3)
    # counts = cxd.get_counts(np.asarray(obs['proton_flux_fit_energy'][idx])*1e3,
    #                         np.asarray(obs['proton_flux_fit'][idx])/1e3,
    #                         svn=72, low=5e3, high=800e3)
    # Integral flux
    omn_gt = interpolate.BSpline(*(interpolate.splrep(fluxmap['energies'], omni_inst)
                                 )).integrate(np.min(fluxlim*1000),np.max(800000))
    fig, axes = plot_omni(cxd, omni_inst, fluxmap_inst, label='ns{} predicted'.format(satnum))
    omni_source = cxd.calculate_omni(fluxmap_source, source=obsdict, initialE=True, from_look=False)
    print(">{}MeV: obs={}; pred={}".format(fluxlim, obs_gt, omn_gt))
    add_extra_omni(axes[0], omni_source, fluxmap_source, label='Source (@ {})'.format(obsname), color='green')
    ## Plot initial spectrum
    # axes[0].plot(obsdict['energy']/1000, obsdict['flux']*1000, 'k-.')
    axes[0].legend()
    ylims = axes[0].get_ylim()
    cdict = cxd.cutoff_info
    axes[0].plot(cdict['ec_low'], ylims[0], marker='^', mec='k', mfc='silver', clip_on=False)
    axes[0].plot(cdict['ec_high'], ylims[0], marker='^', mec='k', mfc='grey', clip_on=False)
    axes[0].set_ylim(ylims)
    plt.savefig('flux_spectrum_ns{}_050500.png'.format(satnum))
    # Plot expected, everything is omni flux _per steradian_
    axes[0].plot(obs['proton_flux_fit_energy'][idx], obs['proton_flux_fit'][idx], 'y:', label='ns{}'.format(satnum))
    plt.savefig('flux_spectrum_ns{}_050500_plus_obs.png'.format(satnum))
    # Now plot directional flux at 50 MeV
    plot_directional_flux(50e3, fluxmap, fluxmap_inst)
    plt.savefig('directional_flux_50MeV_ns{}_050500.png'.format(satnum))
    # Now plot observedd vs predicted counts
    counts = np.array([cc for k, cc in counts.items() if k >= 12])
    plot_counts_pred_obs(counts, obs['rate_proton_measured'][idx])
    plt.savefig('counts_ns{}_comp_{}_050500.png'.format(satnum, obsname))
    plt.show()
