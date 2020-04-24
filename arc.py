"""
Classes for orbital arc (angles-only) observations
"""

import sys
import numpy as np
from datetime import (
    datetime,
    timedelta
)
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from st.tle import ST
from st.site import Site
from algorithm.gauss import gaussAlgorithm
from algorithm.utils import (
    Log,
    simulateEphem
)

class Arc:
    """
    Orbital arc
    """
    def __init__(self, site, ephem_table=None):
        """
        Initialise Arc

        Parameters
        ----------
        site : st.site.Site
            Observation site
        ephem_table : astropy.table.Table, optional
            Table containing 'UTC', 'RA', 'DEC', 'RAERR', 'DECERR' information for Arc
            Default = None (to simulate ephemeris table)
        """
        if not self._validate(site,
                              ephem_table):
            sys.exit()

    def _validate(self, site, ephem_table):
        """
        Validate input

        Parameters
        ----------
        site : st.site.Site
            Observation site
        ephem_table : astropy.table.Table | NoneType
            Table containing 'UTC', 'RA', 'DEC', 'RAERR', 'DECERR' information for Arc | NoneType
        """
        self.site = site
        if not isinstance(self.site, Site):
            print('ArcError: invalid observation site')
            return 0
        self.ephem = ephem_table
        if self.ephem is not None:
            if not isinstance(self.ephem, Table):
                print('ArcError: invalid ephemeris table')
                return 0
            if not ['UTC', 'RA', 'DEC', 'RAERR', 'DECERR'] <= self.ephem.colnames:
                print('ArcError: ephemeris table must contain "UTC", "RA", "DEC", "RAERR", "DECERR"')
                return 0
            self.ephem = self.ephem[np.argsort(self.ephem['UTC'])]
        else:
            print('ArcWarning: no ephemeris table provided')
        return 1

    def set_ephem(self, ephem_table):
        """
        Set ephemeris table for Arc

        Parameters
        ----------
        ephem_table : astropy.table.Table
            Table containing 'UTC', 'RA' and 'DEC' information for Arc
        """
        self.ephem = ephem_table

    def simulate_ephem(self, norad_id, start, end, step, out_dir=None):
        """
        Simulate ephemeris information for a given NORAD object

        Parameters
        ----------
        norad_id : int
            NORAD ID for object
        start : datetime.datetime
            Start epoch
        end : datetime.datetime
            End epoch
        step : datetime.timedelta
            Time step between simulated observations
        out_dir : str, optional
            Path to output directory
            Default = None

        Returns
        -------
        ephem : astropy.table.Table
            Table containing simulated ephemeris information
        """
        return simulateEphem(self.site,
                             norad_id,
                             start,
                             end,
                             step,
                             out_dir=out_dir)

    def initialOrbit(self, idx1, idx2, idx3, improve=False, verbose=False):
        """
        Carry out Gauss algorithm of preliminary orbit determination
        NB: currently assumes central body is Earth

        Parameters
        ----------
        idx1, idx2, idx3 : int
            Ephemeris indices for Gauss method
        improve : bool, optional
            Toggle to carry out algorithm to improve initial orbit
            Default = False
        verbose : bool, optional
            Toggle to display orbital information
            Default = False

        Returns
        -------
        output : dict | None
            Output information with keys:
            tau_1, tau_3 - time intervals
            [1. i] for i orbital solutions, each containing:
                a     - semimajor axis [km]
                r_p   - radius at perihelion [km]
                r_a   - radius at aphelion [km]
                e     - eccentricity
                i     - inclination [degrees]
                Omega - right ascension of the ascending node [degrees]
                omega - argument of perigee [degrees]
                theta - true anomaly [degrees]
                T     - period [hours]
            | None upon failure
        """
        self._set_obs(idx1,
                      idx2,
                      idx3)
        return gaussAlgorithm(self,
                              improve=improve,
                              verbose=verbose)

    def _set_obs(self, idx1, idx2, idx3):
        """
        Set trio of observations for Gauss method

        Parameters
        ----------
        idx1, idx2, idx3 : int
            Ephemeris indices for Gauss method
        """
        self.obs1 = Observation(self.site,
                                self.ephem[idx1])
        self.obs2 = Observation(self.site,
                                self.ephem[idx2])
        self.obs3 = Observation(self.site,
                                self.ephem[idx3])

class Observation:
    """
    Angles-only observation at a particular epoch
    """
    def __init__(self, site, ephem_entry):
        """
        Initialise Observation

        Parameters
        ----------
        site : site.Site
            Observation site
        ephem_entry : astropy.table.Table
            Entry from ephemeris table
        """
        self.site = site
        self.coord = SkyCoord(ra=ephem_entry['RA'],
                              dec=ephem_entry['DEC'],
                              unit=(u.deg, u.deg))
        self.utc = datetime.strptime(ephem_entry['UTC'],
                                     '%Y-%m-%dT%H:%M:%S.%f')
        self._get_ha()

    def _get_lst(self):
        """
        Obtain Local Sidereal Time at epoch of Observation
        """
        self.lst = Time(self.utc,
                        scale='utc',
                        location=self.site.geodetic).sidereal_time('apparent')

    def _get_ha(self):
        """
        Obtain hour angle for Observation
        """
        self._get_lst()
        self.hourangle = (self.lst - self.coord.ra).wrap_at(12 * u.hourangle)

