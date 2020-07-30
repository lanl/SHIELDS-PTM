import os
import glob
import numpy as np
from scipy import constants
from scipy import special
from scipy import linalg
from scipy import integrate
from scipy import interpolate
# local package import
# local package imports
try:
    from . import ptm_tools as pt
except ImportError:
    # if running as a script in source directory
    # relative import won't work
    import ptm_tools as pt


class ptm_postprocessor(object):

    __ckm = constants.speed_of_light/1e3
    __csq = (constants.speed_of_light/1e3)**2
    __elec_restmass = 511.0
    __ccm = constants.speed_of_light*1e2

    __supported_distributions = ('kappa', 'kaprel', 'maxwell', 'juttner')

    # Default parameters for available distribution types, assumes electron species
    __params_kappa = {'density': 1.0e-1, 'energy': 5, 'mass': 1.0, 'kappa': 2.5}
    __params_kaprel = {'density': 1.0e-2, 'energy': 100.0, 'mass': 1.0, 'kappa': 5.0}
    __params_maxwell = {'density': 1.0e-1, 'energy': 10, 'mass': 1.0}
    __params_juttner = {'density': 4.73e-3, 'energy': 160.0, 'mass': 1.0}

    def __init__(self, source_type='kappa', params={}, filedir=None):
        """
        """
        if filedir is not None:
            self.__filedir = os.path.expanduser(filedir)
        else:
            self.__filedir = os.getcwd()

        # Initialize postprocessor. Default is kappa distribution but this can
        # be changed using the constructor call or with an explicit call to the
        # configure_source method after the postprocessor is initialized
        self.set_source(source_type, params)


    def set_source(self,source_type='kappa', params={}):
        """Set source distribution type and parameters

        Optional Arguments
        ------------------
        source_type : str
            Distribution name. Default 'kappa', other options are
            'kaprel', 'maxwell', 'juttner'
        params : dict
            Distribution parameters, all distributions have 'density',
            'energy', and 'mass'. The two kappa distributions also have
            'kappa'. Mass defaults to 1 in all cases, and is given in
            units of the electron mass.

        Distribution Parameters
        -----------------------
        kappa : (non-relativistic) kappa distribution
            density - number density, default 1e-1 cm^-2
            energy - characteristic energy, default 5 keV
            kappa - kappa parameter, default 2.5
        kaprel : relativistic Kappa-type distribution
            density - number density, default 1e-2 cm^-2
            energy - characteristic energy, default 100 keV
            kappa - kappa parameter, default 5.0
        maxwell : Maxwell-Boltzmann distribution
            density - number density, default 1e-1
            energy - characteristic energy, default 10 keV
        juttner : Maxwell-Juttner (relativistic) distribution
            density - number density, default 4.73e-3
            energy - characteristic energy, default 160 keV

        Notes
        -----
        The asymptotic characteristic energies differ between 'kappa' and 'kaprel',
        so for comparison in the low-energy regime Ec(kappa) ~ 2*Ec(kaprel).
        
        """
        if source_type not in self.__supported_distributions:
            raise ValueError('Specified source type ({:}) not supported'.format(source_type))
        elif source_type == 'kappa':
            self.__use_kappa(params)
            self.get_flux = self.__flux_kappa
            self.get_dist = self.__dist_kappa
            self.get_dist_u = self.__dist_kappa_u
        elif source_type == 'kaprel':
            self.__use_kaprel(params)
            self.get_flux = self.__flux_kaprel
            self.get_dist = self.__dist_kaprel
            self.get_dist_u = self.__dist_kaprel_u
        elif source_type == 'maxwell':
            self.__use_maxwell(params)
            self.get_flux = self.__flux_maxwell
            self.get_dist = self.__dist_maxwell
            self.get_dist_u = self.__dist_maxwell_u
        elif source_type == 'juttner':
            self.__use_juttner(params)
            self.get_flux = self.__flux_juttner
            self.get_dist = self.__dist_juttner
            self.get_dist_u = self.__dist_juttner_u

        self.source_type = source_type


    def __use_kappa(self, params={}):
        # Set parameters of the postprocessor for non-relativistic kappa option
        # Get the parameters by merging user specifications with defaults
        my_params = {**self.__params_kappa, **params}

        if my_params['kappa'] <= 1.5:
            raise ValueError('__use_kappa: kappa < 1.5 ({:}) is not valid'.format(my_params['kappa']))

        self.__kappa = my_params['kappa']
        self.__n = my_params['density']
        self.__ec = my_params['energy']
        self.__mc2 = my_params['mass']*self.__elec_restmass
        self.__wc = self.__ckm*np.sqrt((2*self.__kappa-3)*self.__ec/(self.__kappa*self.__mc2))
        self.__fcoef = self.__n*self.__wc*(self.__ccm/self.__ckm)*np.sqrt(4*np.pi)/(self.__ec*self.__kappa*special.beta(self.__kappa-1/2,3/2))
        self.__dcoef = self.__n/(2*np.pi*(self.__kappa*self.__wc**2)**(3/2)*special.beta(self.__kappa-1/2,3/2))


    def __use_kaprel(self, params={}):
        # Set parameters of the postprocessor for relativistic kappa option
        # Get the parameters by merging user specifications with defaults
        my_params = {**self.__params_kappa, **params}

        if my_params['kappa'] <= 1.5:
            raise ValueError('__use_kaprel: kappa < 1.5 ({:}) is not valid'.format(my_params['kappa']))

        self.__kappa = my_params['kappa']
        self.__n = my_params['density']
        self.__ec = my_params['energy']
        self.__mc2 = my_params['mass']*self.__elec_restmass
        self.__energy_ratio = self.__mc2/self.__ec

        a = 4*np.pi
        b = 8*special.beta(3/2,self.__kappa-2)/(2*self.__kappa-1)
        x = 1-2*self.__energy_ratio/self.__kappa
        c = special.hyp2f1(self.__kappa+1,5/2, self.__kappa+1/2,x)*3
        d = special.hyp2f1(self.__kappa+1,3/2, self.__kappa+1/2,x)*(self.__kappa-2)

        self.__fcoef = self.__n*self.__ccm/(self.__ec*a*b*(c+d))
        self.__dcoef = self.__n/(a*b*(c+d)*self.__ckm**3)


    def __use_maxwell(self, params={}):
        # Set parameters of the postprocessor for non-relativistic Maxwellian option
        # Get the parameters by merging user specifications with defaults
        my_params = {**self.__params_juttner, **params}

        self.__n = my_params['density']
        self.__mc2 = my_params['mass']*self.__elec_restmass
        self.__ec = my_params['energy']
        self.__energy_ratio = self.__mc2/self.__ec
        self.__dcoef = self.__n*(self.__energy_ratio/(2*np.pi*self.__ckm**2))**(3/2)
        self.__fcoef = self.__n*self.__ccm*(self.__energy_ratio/(2*np.pi))**(3/2)/self.__mc2


    def __use_juttner(self, params={}):
        # Set parameters of the postprocessor for Maxwell-Juttner option
        # Get the parameters by merging user specifications with defaults
        my_params = {**self.__params_juttner, **params}

        self.__n = my_params['density']
        self.__mc2 = my_params['mass']*self.__elec_restmass
        self.__ec = my_params['energy']
        self.__energy_ratio = self.__mc2/self.__ec
        self.__dcoef = self.__n*self.__energy_ratio/(4*np.pi*self.__ckm**3*special.kve(2,self.__energy_ratio))
        self.__fcoef = self.__n*self.__ccm*self.__energy_ratio/(4*np.pi*self.__mc2*special.kve(2,self.__energy_ratio))


    def __dist_kappa(self, energy):
        # Kappa distribution function that takes energy as input
        v = np.sqrt(2*energy/self.__mc2)
        fun = 1/(1+(v/(self.__wc/self.__ckm))**2/self.__kappa)**(self.__kappa+1)
        f = self.__dcoef*fun

        return f


    def __dist_kappa_u(self, v):
        # Kappa distribution function that takes velocity as an input
        # v=p/mc is the normalized momentum, presumed non-relativistic
        fun = 1/(1+(v/(self.__wc/self.__ckm))**2/self.__kappa)**(self.__kappa+1)
        f = self.__dcoef*fun

        return f


    def __dist_kaprel(self, energy):
        # Relativistic kappa distribution that takes energy as an input
        gam = 1+energy/self.__mc2
        fun = 1/(1+(gam-1)*self.__energy_ratio/self.__kappa)**(self.__kappa+1)
        f = self.__dcoef*fun

        return f


    def __dist_kaprel_u(self, u):
        # Relativistic kappa distribution that takes momentum as an input
        # u=p/mc is the normalized momentum
        gam = np.sqrt(1+u**2)
        fun = 1/(1+(gam-1)*self.__energy_ratio/self.__kappa)**(self.__kappa+1)
        f = self.__dcoef*fun

        return f


    def __dist_maxwell(self, energy):
        # Non-relativistic Maxwellian distribution
        beta2 = 2*energy/self.__mc2
        fun = np.exp(-self.__energy_ratio*beta2/2)
        f = self.__dcoef*fun

        return f


    def __dist_maxwell_u(self, u):
        # Non-relativistic Maxwellian distribution that takes momentum as an input
        # u = p/mc = v/c in non-relativistic limit
        beta2 = u*u
        fun = np.exp(-self.__energy_ratio*beta2/2)
        f = self.__dcoef*fun

        return f


    def __dist_juttner(self, energy):
        # Maxwell-Juttner distribution that takes energy as an input
        gam = 1 + energy/self.__mc2
        fun = np.exp(-(gam-1)*self.__energy_ratio)
        f = self.__dcoef*fun

        return f


    def __dist_juttner_u(self, u):
        # Maxwell-Juttner distribution that takes momentum as an input
        # u=p/mc is the normalized momentum
        gam = np.sqrt(1+u**2)
        fun = np.exp(-self.__energy_ratio*(gam-1))
        f = self.__dcoef*fun

        return f


    def __flux_kappa(self, energy):
        # Based on equation 10a of Vasyliunas [1962]. The coefficient is caluclated in
        # the __use_kappa method
        x = energy/self.__ec
        fun = x/(1+x/self.__kappa)**(self.__kappa+1)

        j = self.__fcoef*fun

        return j


    def __flux_kaprel(self, energy):
        # Based on Equations 4-5 of Xiao, Zhou, Li, and Cai [2008]. The coefficent
        # is calculated in the __use_kaprel method
        gamma = 1+energy/self.__mc2
        fun = (gamma*gamma-1)/(1+(gamma-1)*self.__energy_ratio/self.__kappa)**(self.__kappa+1)
        j = self.__fcoef*fun

        return j


    def __flux_maxwell(self, energy):
        # Standard Maxwell-Boltzmann distribution
        beta2 = 2*energy/self.__mc2
        fun = beta2*np.exp(-self.__energy_ratio*beta2/2)

        j = self.__fcoef*fun

        return j


    def __flux_juttner(self, energy):
        # Based on Equation 14 of Morley, Sullivan, Schirato, and Terry [2014]. The coefficent
        # is calculated in the __use_juttner method
        gamma = 1 + energy/self.__mc2
        fun = (gamma*gamma-1)*np.exp(-(gamma-1)*self.__energy_ratio)

        j = self.__fcoef*fun

        return j


    def map_flux(self, fluxmap):
        """
        Map fluxes from source region to observation point using Liouville's theorem
        """
        # Energies
        ei = fluxmap['init_E']
        ef = fluxmap['final_E']

        # Calculate the flux at the source
        ji = self.get_flux(ei)

        # Relativistic gammas
        gami = 1 + ei/self.__mc2
        gamf = 1 + ef/self.__mc2

        # Mapped the flux to the observation point, the prefactor is ratio of squared momenta
        jf = ((gamf*gamf-1)/(gami*gami-1))*ji

        return jf


    def get_omni_flux(self, pa, flux, method='spline'):
        """
        This features two methods for omnidirectional flux calculations. The default is to fit B-splines to the
        fluxes and then integrate the spline functions. Alternatively, by passing "method = bin", a trapezoidal
        integration scheme can be used.
        """
        fshape = np.shape(flux)
        npa = np.size(pa)

        # Infer if the pitch angles are radians or degrees
        if np.max(np.abs(pa)) > np.pi:
            mya = np.radians(pa)
        else:
            mya = pa

        # Determine which dimension corresponds to pitch angles
        if fshape[1] == npa:
            myflux = flux
        elif fshape[0] == npa:
            myflux = np.transpose(flux)
        else:
            raise ValueError("Dimension mismatch in get_omni_flux")

        if method == 'spline':
            # Calculate the omnidiretional flux using B-spline integration, with 2*pi accounting for gyrophase
            sina = np.sin(mya)
            omni = 2*np.pi*np.array([interpolate.BSpline(*(interpolate.splrep(mya,f*sina))).integrate(np.min(mya),np.max(mya)) for f in myflux])
        else:
            coef = 2*np.pi*(mya[1:] - mya[:-1])
            favg = 0.5*(myflux[:, 1:] + myflux[:, :-1])
            savg = np.sin(0.5*(mya[1:] + mya[:-1]))

            # This reduction approximates the weighted integral over pitch angles
            omni = np.einsum("i,ji", coef*savg, favg)

        return omni


    def process_run(self, runid, verbose=True):
        """
        Read in the results of a PTM simulation and calculate the fluxes

        -------
        Inputs
        -------

        runid : integer, list, or None
            Integer - identification tag for run to be analyzed (e.g. 1 for data in map_0001.dat)
            List - List of integer ID tags for runs to be combined and analyzed
            None - Use all map_*.dat files in the selected directory

        -------
        Outputs
        -------

        results     dictionary  fluxes and associated quantities describing the ptm fluxes. Has keys:
                                "fluxmap"   Raw fluxmap from file
                                "energies"  Energies at which fluxes were calculated
                                "angles"    Pitch angles at which fluxes were calculated
                                "flux"      Differential fluxes
                                "omni"      Omnidirectional fluxes
                                "kappa"     Spectral index of kappa distribution
                                "n_dens"    Particle density in source region
                                "e_char"    Characteristic energy of kappa distribution

        """
        results = {}

        if runid is None:
            fname = glob.glob(os.path.join(self.__filedir, 'map_*.dat'))
        else:
            try:
                # handle an iterable of run IDs
                fname = [os.path.join(self.__filedir, 'map_{0:04}.dat'.format(ri))
                         for ri in runid]
            except IndexError:
                # fall back to single integer provided
                fname = self.__filedir+'/map_{:04}.dat'.format(runid)

        fluxmap = pt.parse_map_file(fname)

        flux = self.map_flux(fluxmap)
        omni = self.get_omni_flux(fluxmap['angles'],flux)

        results['position'] = fluxmap['final_x']
        # results['fluxmap'] = fluxmap
        results['initial_E'] = fluxmap['init_E']
        results['final_E'] = fluxmap['final_E']
        results['energies'] = fluxmap['energies']
        results['angles'] = fluxmap['angles']
        results['flux'] = flux
        results['omni'] = omni

        results['kappa'] = self.__kappa
        results['n_dens'] = self.__n
        results['e_char'] = self.__ec

        if verbose:
            print("Energy grid : ", fluxmap['energies'])
            print("PitchAngle grid : ", fluxmap['angles'])
            print("Final Particle Energies [PA] : ",  fluxmap['final_E'])
            print("Diff Flux [E[PA]]: ", flux)
            print("Omni Flux [E]: ", omni)

        return results


    def test_distributions(self):
        """
        Check valididty of distribution functions.
        """
        # First, check that velocity/momentum and energy representations give the same answer.
        # Second, confirm that distribution functions are properly normalized. Each gets multiplied by 4 pi c^3 / n. Note that we are
        # using the versions that take normalized momentum as the argument rather than energy b/c this simplifies the
        # mathematics. The other versions yield the same results when momentum is converted to kinetic energy

        # Store current parameters
        params = {'energy':self.__ec,'density':self.__n,'mass':self.__mc2/self.__elec_restmass,'source':self.source_type}

        # Check if we're currently in the kappa family, and if so remember to store kappa as well
        if hasattr(self,'_ptm_postprocessor__kappa'):
            params['kappa'] = self.__kappa

        c3 = self.__ckm**3

        # Maxwell distribution
        # Integrate to show proper normalization
        self.set_source('maxwell')
        res0=integrate.quad(lambda x:4*np.pi*c3*x*x*self.__dist_maxwell_u(x),0,np.inf)[0]/self.__n

        # Show equivalence between distributions
        test_speed=0.1 #Normalized momentum p/mc
        energy = (1/2)*self.__mc2*(test_speed**2)
        m1,m2 = self.get_dist_u(test_speed), self.get_dist(energy)

        #Maxwell-Juttner distribution
        # Integrate to show proper normalization
        self.set_source('juttner')
        res1=integrate.quad(lambda x:4*np.pi*c3*x*x*self.__dist_juttner_u(x),0,np.inf)[0]/self.__n

        # Show equivalence between distributions
        test_momentum=0.1 #Normalized momentum p/mc
        gamma=np.sqrt(1+test_momentum**2)
        energy=(gamma-1)*self.__mc2 # E = (gamma-1)*mc^2
        j1,j2 = self.get_dist_u(test_momentum), self.get_dist(energy)

        # Integrate Kappa distribution
        self.set_source('kappa')
        res2=integrate.quad(lambda x:4*np.pi*c3*x*x*self.__dist_kappa_u(x),0,np.inf)[0]/self.__n

        # Show equivalence between distributions using energy and velocity inputs
        test_speed=0.1 #Fraction of light speed
        energy=(1/2)*self.__mc2*test_speed**2 # E = (1/2)(mc^2)*(v/c)**2
        k1,k2 = self.get_dist_u(test_speed),self.get_dist(energy)

        # Integrate relativistic Kappa distribution
        self.set_source('kaprel')
        res3=integrate.quad(lambda x:4*np.pi*c3*x*x*self.__dist_kaprel_u(x),0,np.inf)[0]/self.__n

        gamma=np.sqrt(1+test_momentum**2)
        energy=(gamma-1)*self.__mc2 # E = (gamma-1)*mc^2

        r1,r2 = self.get_dist_u(test_momentum),self.get_dist(energy)

        print('USING DEFAULT ARGUMENTS\n')

        print('Maxwell Distribution:')
        print('Compare f(u) and f(E) (should be same):')
        print('{:}\t{:}'.format(m1,m2))
        print('Integral of distribution function (should be 1):')
        print('{:}\n'.format(res0))

        print('Juttner Distribution:')
        print('Compare f(u) and f(E) (should be same):')
        print('{:}\t{:}'.format(j1,j2))
        print('Integral of distribution function (should be 1):')
        print('{:}\n'.format(res1))

        print('Kappa Distribution:')
        print('Compare f(v) and f(E) (should be same):')
        print('{:}\t{:}'.format(k1,k2))
        print('Integral of distribution function (should be 1):')
        print('{:}\n'.format(res2))

        print('Relativisitc Kappa Distribution:')
        print('Compare f(u) and f(E) (should be same):')
        print('{:}\t{:}'.format(r1,r2))
        print('Integral of distribution function (should be 1):')
        print('{:}\n'.format(res3))

        # Return to the original source type
        self.set_source(params['source'], params)


    def test_omni(self):
        """
        Calculate omnidirectional flux for simple test cases
        """
        # Use sin^n(alpha) as the PAD and assume gyrotropic (2*pi factor)
        import matplotlib.pyplot as plt

        def notch(x,n):
            # Notched sine function to approximate missing data point
            res = np.sin(x)**n*(1-np.exp(-((x-np.pi/2)/0.01)**2))
            return res

        omni = np.array([2*np.pi*integrate.quad(lambda x: np.sin(x)*notch(x,n), 0, np.pi)[0]
                         for n in range(1,6)])

        # First, uniformly sampled in angle space
        qtest = np.linspace(0, np.pi, 13)
        ftest = np.array([notch(qtest, n) for n in range(1, 6)])
        omnia = self.get_omni_flux(qtest, ftest)
        plt.plot(qtest, ftest[0, :], 'ro', ms=10)
        plt.plot(qtest, ftest[1:, :].T, 'ro', ms=10, label='_nolegend_')

        # Next, uniformly sampled in cos(pa) space
        qtest = np.arccos(np.linspace(1, -1, 13))
        ftest = np.array([notch(qtest, n) for n in range(1, 6)])
        omnic = self.get_omni_flux(qtest, ftest)
        plt.plot(qtest, ftest[0, :], 'ko', ms=6)
        plt.plot(qtest, ftest[1:, :].T, 'ko', ms=6, label='_nolegend_')

        print("Omnidirection flux calculations:\n")

        print("Using spline method:")
        print('{:}\t{:}'.format('Full Integral', omni))
        print('{:}\t{:}'.format('Uniform PA', omnia))
        print('{:}\t{:}\n'.format('Uniform COS(PA)', omnic))

        # First, uniformly sampled in angle space
        qtest = np.linspace(0, np.pi, 13)
        ftest = np.array([notch(qtest, n) for n in range(1, 6)])
        omnia = self.get_omni_flux(qtest, ftest, method='bin')

        # Next, uniformly sampled in cos(pa) space
        qtest = np.arccos(np.linspace(1, -1, 13))
        ftest = np.array([notch(qtest, n) for n in range(1, 6)])
        omnic = self.get_omni_flux(qtest, ftest, method='bin')

        print("Using binning method:")
        print('{:}\t{:}'.format('Full Integral', omni))
        print('{:}\t{:}'.format('Uniform PA', omnia))
        print('{:}\t{:}\n'.format('Uniform COS(PA)', omnic))

        xv = np.linspace(0, np.pi, 1000)
        for n in range(1, 6):
            yv = notch(xv, n)
            plt.plot(xv, yv)

        plt.title(r'$f(\alpha)=(1-\exp(-100(\alpha-\pi/2)^2)\sin^{n}(\alpha)$', fontsize=16)
        plt.legend(['Uniform PA','Uniform cos(PA)','n=1','n=2','n=3','n=4','n=5'],
                   frameon=False, loc='upper left')


def process_ram_boundary(self, griddir=None, write_files=True, outdir=None,
                         date={'year': 2000, 'month': 1, 'day': 1}):
    """
    Rungrids are a new configuration option provided in the updated version of ptm_input.
    In addition to the input files, there is a rungrid.txt file that
    describes the characteristics of each input file (time and location).

    When the RAM boundary is simulated using PTM via rungrid configuration, we are able to
    simplify the post-processing workflow using this routine.
    """

    if griddir == None:
        mydir = self.__filedir
    else:
        mydir = griddir

    fname = mydir+'/rungrid.txt'

    if os.path.isfile(fname):
        rungrid = np.loadtxt(fname, skiprows=1)

        # This is backwards tracing, so fluxes will be calculated at the later
        # time, which is given in the third column

        runids = map(int,rungrid[:,0])
        times = np.unique(rungrid[:,2])
        rvals = np.unique(rungrid[:,3])

        if not np.allclose(rvals, rvals[0], 1e-3):
            raise Exception('Error in process_ram_boundary: points are not at fixed radial distance.')

        fluxdata = {'rungrid': rungrid}
        fluxdata['runid'] = runids
        fluxdata['times'] = times
        fluxdata['R'] = rvals[0]
        fluxdata['mlt'] = np.sort(np.unique(rungrid[:,4]))

        for runid in runids:
            fluxdata[runid] = self.process_run(runid)
    else:
        raise Exception('Error in process_rungrid: ' + fname + ' not found')

    if(write_files):
        self.write_ram_fluxes(fluxdata, date=date, outdir=outdir)

    return fluxdata


def write_ram_fluxes(self, fluxdata, date={'year': 2000, 'month': 1, 'day': 1}, outdir=None):
    """
    Write time- and space-dependent fluxes in a RAM boundary file
    """
    cadence = (fluxdata['times'][1]-fluxdata['times'][0])//60

    year = date['year']
    month = date['month']
    day = date['day']

    fname = '{:4}{:02}{:02}_ptm_geomlt_{:}-min.txt'.format(year, month, day, cadence)

    nenergy = fluxdata[1]['energies'].size

    # Header lines
    header1 = '# PTM Particle Fluxes for RAM\n'
    header2 = '# Header Format string: (a24,a6,2x,a72,36a18)\n'
    header3 = '# DATA   Format string: (a24,f6.1,2x,36(i2),36(f18.4))\n'

    # This parameter has to be in the file but it's not used by RAM
    nsc = np.ones([nenergy],dtype='int')
    nscstring = (nenergy*'{:2}').format(*nsc)

    # Formatting strings
    headFormat = '{:>24}{:>6}  {:>72}'
    dataFormat = '{:6.1f}  ' + nscstring + (nenergy*'{:18.4f}')
    timeFormat = '{:4}-{:02}-{:02}T{:02}:{:02}:{:02}.000Z'

    with open(fname,'w') as fh:
        fh.writelines(header1);
        fh.writelines(header2);
        fh.writelines(header3);
        fh.writelines(headFormat.format('CCSDS','MLT','NSC')+(nenergy*'{:18}'+'\n').format(*fluxdata[1]['energies']))

        i=0
        for time in fluxdata['times']:
            hour,minute,second = self.seconds_to_hhmmss(time)
            for mlt in fluxdata['mlt']:
                i += 1
                saflux = fluxdata[i]['omni']/(4.0*np.pi)
                # Note asterisk in format statement (*np.r_), this is required for correct passing of values
                dataline = timeFormat.format(year, month, day, hour, minute, second)\
                           + dataFormat.format(*np.r_[mlt, saflux]) + '\n'
                fh.writelines(dataline)


def tm03_moments(x, y, sw_user=None):
    """
    The Tsyganeko and Mukai [2003] plasma sheet model (with corrected equations).

    Inputs:
      x     GSM x position in Earth radii
      y     GSM y position in Earth radii
      sw_user   An optional dictionary which may contain the following keys:
          'bperp' Perpendicular component of the solar wind magnetic field in nT
          'theta' Solar wind clock angle in degrees
          'vx'    Solar wind speed in km/s
          'n'     Solar wind density in cm-3
          'p'     Solar wind dynamic pressure in nPa

    Default values of swD are {'bperp':5,'theta':90,'vx':500,'n':10,'p':3}, and the
    user needs only to specify values that are different from these.

    Outputs:
      res   A dictionary with the following keys:
          'P'   Local pressure in nPa
          'T'   Local temperature in keV
          'n'   Local ion/electron density in cm-3

    Reference:

    Tsyganeko, N. A. and T. Mukai (2003), Tail plasma sheet models derived from Geotail particle data,
        J. Geophys. Res., 108(A3),1136, doi:10.1029/2002JA009707

    The output from the TM03 model can be used to set n and Ec inputs to the energy_to_flux routine.

    Jesse Woodroffe
    6/2/2016

    Modified to simplify interface and eliminate "double-call" initialization
    """
    TM03_defaults = {'bperp': 5.0, 'theta': 90.0, 'vx': 500.0, 'n': 10.0, 'p': 3.0}

    swD = {**TM03_defaults} if sw_user is None else {**TM03_defaults,**sw_user}

    # Evaluate TM03 model using provided parameters
    aT = np.r_[ 0.0000, 1.6780,-0.1606, 1.6690, 4.8200, 2.8550,-0.6020,-0.8360,
               -2.4910, 0.2568, 0.2249, 0.1887,-0.4458,-0.0331,-0.0241,-2.6890,
                1.2220]
    aN = np.r_[ 0.0000,-0.1590, 0.6080, 0.5055, 0.0796, 0.2746, 0.0361,-0.0342,
               -0.7935, 1.1620, 0.4756, 0.7117]
    aP = np.r_[ 0.0000, 0.0570, 0.5240, 0.0908, 0.5270, 0.0780,-4.4220,-1.5330,
               -1.2170, 2.5400, 0.3200, 0.7540, 1.0480,-0.0740, 1.0150]

    bperp = swD['bperp']/5
    bz = bperp*np.cos(np.deg2rad(swD['theta']))

    bzn, bzs = (bz, 0) if bz > 0 else (0, -bz)

    vsw = swD['vx']/500.0
    nsw = swD['n']/10.0
    fsw = bperp*np.sqrt(np.sin(np.deg2rad(swD['theta'])/2))
    rho = np.sqrt(x*x + y*y)/10.0
    psw = swD['p']/3.0
    phi = -np.arctan2(y, x)
    rm1 = rho-1.0

    T = (aT[1]*vsw + aT[2]*bzn+aT[3]*bzs+aT[4]
         * np.exp(-(aT[9]*vsw**aT[15] + aT[10]*bzn+aT[11]*bzs)*rm1)
         + (aT[5]*vsw+aT[6]*bzn+aT[7]*bzs+aT[8]*
         np.exp(-(aT[12]*vsw**aT[16]+aT[13]*bzn+aT[14]*bzs)*rm1))*np.sin(phi)**2)

    N = ((aN[1]+aN[2]*nsw**aN[10]+aN[3]*bzn+aN[4]*vsw*bzs)*rho**aN[8]+
         (aN[5]*nsw**aN[11]+aN[6]*bzn+aN[7]*vsw*bzs)*rho**aN[9]*np.sin(phi)**2)

    P = (aP[1]*rho**aP[6] + aP[2]*psw**aP[11]*rho**aP[7]
         + aP[3]*fsw**aP[12]*rho**aP[8] + (aP[4]*psw**aP[13]*np.exp(-aP[9]*rho)
         + aP[5]*fsw**aP[14]*np.exp(aP[10]*rho)) *np.sin(phi)**2)

    res = {'P': P, 'T': T, 'n': N}

    return res


if __name__ == "__main__":
    pp = ptm_postprocessor()
    pp.test_distributions()
    pp.test_omni()
