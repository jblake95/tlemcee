"""
Config reading for tlemcee
"""

import os
import sys
import json
import numpy as np
from datetime import datetime
from astropy.table import Table
from skyfield.api import utc

from arc import Arc
from st.tle import TLE
from st.site import Site

class Config:
    """
    Container for tlemcee config
    """
    def __init__(self, config_path):
        """
        Initialise Config

        Parameters
        ----------
        config_path : str
            Path to config JSON file
        """
        if self._validate(config_path):
            self._read_config()
        else:
            sys.exit()

    def _validate(self, config_path):
        """
        Validate input

        Parameters
        ----------
        config_path : str
            Path to config JSON file
        """
        self._input = config_path
        if not os.path.exists(self._input):
            print('ConfigError: Input path not found')
            return 0
        return 1

    def _read_config(self):
        """
        Read config information
        """
        with open(self._input, 'r') as cf:
            self._config = json.load(cf)

        # data information
        self._init_data()

        # prior information
        self._init_priors()

        # sampler information
        self._init_sampler()

    def _init_data(self):
        """
        Read in data from ephemeris file
        """
        if 'site_name' in self._config.keys():
            self.site = Site(self._config['site_name'])
        else:
            print('ConfigWarning: No site_name specified, defaulting to LaPalma')
            self.site = Site('LaPalma')

        if 'ephem_path' in self._config.keys():
            self.ephem_path = self._config['ephem_path']
            if os.path.exists(self.ephem_path):
                arc = Arc(self.site, Table.read(self.ephem_path, format='csv'))
            else:
                print('ConfigError: Input ephemeris path not found')
                sys.exit()
        else:
            print('ConfigError: No ephemeris file specified')
            sys.exit()

        times = []
        for time in arc.ephem['UTC']:
            times.append(datetime.strptime(time, '%Y-%m-%dT%H:%M:%S.%f'))

        self.data = {'time': np.array(times),
                     'ra': np.array(arc.ephem['RA']),
                     'dec': np.array(arc.ephem['DEC']),
                     'raerr': np.array(arc.ephem['RAERR']),
                     'decerr': np.array(arc.ephem['DECERR'])}

        # pre-calculate repeatedly used terms in lnlike
        self.data.update({'ra_sigma2': self.data['raerr'] ** 2})
        self.data.update({'dec_sigma2': self.data['decerr'] ** 2})
        self.data.update({'ra_lnsigma2': np.log(self.data['ra_sigma2'])})
        self.data.update({'dec_lnsigma2': np.log(self.data['dec_sigma2'])})

        if 'out_dir' in self._config.keys():
            self.out_dir = self._config['out_dir']
            if not os.path.exists(self.out_dir):
                try:
                    os.mkdir(self.out_dir)
                except:
                    print('ConfigWarning: Failed to set up output directory')
                    self.out_dir = None
        else:
            print('ConfigWarning: No output directory specified')
            self.out_dir = None

    def _init_tle(self):
        """
        Construct input TLE from fixed priors
        """
        self.tle = TLE()  # blank canvas

        # set fixed priors if any
        if 'fixed' in self._config.keys():
            count = 0
            fixed_params = self._config['fixed']
            if 'epoch' in fixed_params:
                epoch = datetime.strptime(fixed_params['epoch']['value'],
                                          '%Y-%m-%dT%H:%M:%S.%f')
                count += 1
            else:
                epoch = None
            if 'mmdot' in fixed_params:
                mmdot = float(fixed_params['mmdot']['value'])
                count += 1
            else:
                mmdot = None
            if 'mmdot2' in fixed_params:
                mmdot2 = float(fixed_params['mmdot2']['value'])
                count += 1
            else:
                mmdot2 = None
            if 'drag' in fixed_params:
                drag = float(fixed_params['drag']['value'])
                count += 1
            else:
                drag = None
            if 'inclination' in fixed_params:
                inclination = float(fixed_params['inclination']['value'])
                count += 1
            else:
                inclination = None
            if 'raan' in fixed_params:
                raan = float(fixed_params['raan']['value'])
                count += 1
            else:
                raan = None
            if 'eccentricity' in fixed_params:
                eccentricity = float(fixed_params['eccentricity']['value'])
                count += 1
            else:
                eccentricity = None
            if 'argperigee' in fixed_params:
                argperigee = float(fixed_params['argperigee']['value'])
                count += 1
            else:
                argperigee = None
            if 'meananomaly' in fixed_params:
                meananomaly = float(fixed_params['meananomaly']['value'])
                count += 1
            else:
                meananomaly = None
            if 'mm' in fixed_params:
                mm = float(fixed_params['mm']['value'])
                count += 1
            else:
                mm = None

            if count != len(fixed_params):
                print('ConfigWarning: Ignoring unrecognised param in fixed priors: {}'.format(param))

            self.tle.modify_elements(epoch=epoch,
                                     mmdot=mmdot,
                                     mmdot2=mmdot2,
                                     drag=drag,
                                     inclination=inclination,
                                     raan=raan,
                                     eccentricity=eccentricity,
                                     argperigee=argperigee,
                                     meananomaly=meananomaly,
                                     mm=mm)
        else:
            print('ConfigWarning: No fixed priors specified')

    def _init_priors(self):
        """
        Read in prior information
        """
        # fixed priors can be fed into sampler via input TLE
        self._init_tle()

        # uniform priors will be fed into sampler via arrays
        if 'uniform' in self._config.keys():
            self.priors = {'name': [],
                           'init': [],
                           'llim': [],
                           'ulim': [],
                           'wght': []}
            uniform = self._config['uniform']
            for param in ['mmdot', 'mmdot2', 'drag', 'inclination', 'raan',
                          'eccentricity', 'argperigee', 'meananomaly', 'mm']:
                if param in uniform:
                    self.priors['name'].append(param)
                    self.priors['init'].append(float(uniform[param]['init']))
                    self.priors['llim'].append(float(uniform[param]['llim']))
                    self.priors['ulim'].append(float(uniform[param]['ulim']))
                    if self.priors['llim'][-1] > self.priors['ulim'][-1]:
                        print('ConfigError: unphysical limits for {}'.format(param))
                    self.priors['wght'].append(float(uniform[param]['wght']))
                    if self.priors['wght'][-1] > (self.priors['ulim'][-1] - self.priors['llim'][-1]):
                        print('ConfigError: unsuitable weight for {}'.format(param))

            if len(self.priors['init']) != len(uniform):
                print('ConfigWarning: Ignoring unrecognised param(s) in uniform priors')
        else:
            print('ConfigError: No uniform priors specified - please rectify')
            sys.exit()

    def _init_sampler(self):
        """
        Read in sampler information
        """
        if 'n_steps' in self._config.keys():
            self.n_steps = self._config['n_steps']
        else:
            print('ConfigWarning: No n_steps specified, defaulting to 1000')
            self.n_steps = 1000

        # recommended walkers is 4 * n_parameters - scaling factor allows for more
        self.n_dim = len(self.priors['init'])
        if 'walker_scaling' in self._config.keys():
            self.walker_scaling = self._config['walker_scaling']
        else:
            print('ConfigWarning: No walker_scaling specified, defaulting to 1')
            self.walker_scaling = 1
        self.n_walkers = 4 * self.n_dim * self.walker_scaling

        self.init_pos = []
        for i in range(self.n_walkers):
            self.init_pos.append(self.priors['init'] + self.priors['wght'] * np.random.randn(self.n_dim))
