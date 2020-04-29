"""
Two-line element (TLE) manipulation
"""

import sys
import numpy as np
import getpass as gp
from spacetrack import SpaceTrackClient
from skyfield.sgp4lib import EarthSatellite
from skyfield.api import (
    load,
    utc
)
from datetime import timedelta
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

from st.tlelements import (
    unpackElements,
    modifyElements
)
from st.utils import parsePropagationInfo

UN = 'J.Blake@warwick.ac.uk'  # change username here
TS = load.timescale()         # save repeated use in iterative loops

class ST:
    """
    Space-Track Interface
    """
    def __init__(self):
        """
        Connects user to Space-Track account
        """
        self.client = SpaceTrackClient(identity=UN,
                                       password=self._request_access())

    def _request_access(self):
        """
        Obtain access details - requests user's password

        Returns
        -------
        pw : str
            Password entered by user
        """
        return gp.getpass('Space-Track password: ')

    def get_latest_tle(self, norad_id):
        """
        Obtain latest TLE for NORAD object

        Parameters
        ----------
        norad_id : int
            NORAD ID for object

        Returns
        -------
        tle : st.tle.TLE
            Latest TLE for NORAD object
        """
        el_set = [line for line in self.client.tle_latest(norad_cat_id=norad_id,
                                                          iter_lines=True,
                                                          ordinal=1,
                                                          format='3le')]
        return TLE(el_set[1],
                   el_set[2],
                   name=el_set[0])

    def get_past_tle(self, norad_id, epoch, search_radius=10):
        """
        Obtain past TLE (closest to given epoch) for NORAD object

        Parameters
        ----------
        norad_id : int
            NORAD ID for object
        epoch : datetime.datetime
            Desired epoch
        search_radius : int, optional
            Search radius either side of epoch [days]
        """
        epoch_range = '{}--{}'.format((epoch - timedelta(days=search_radius)).strftime('%Y-%m-%d'),
                                      (epoch + timedelta(days=search_radius)).strftime('%Y-%m-%d'))

        el_sets = [line for line in self.client.tle(norad_cat_id=norad_id,
                                                    epoch=epoch_range,
                                                    iter_lines=True,
                                                    format='3le')]
        delta_t = []
        for i in range(int(len(el_sets) / 3)):
            tle = TLE(el_sets[3 * i + 1],
                      el_sets[3 * i + 2],
                      name=el_sets[3 * i])
            delta_t.append(abs(tle.epoch.date - epoch))

        idx = np.argmin(delta_t)
        return TLE(el_sets[3 * idx + 1],
                   el_sets[3 * idx + 2],
                   el_sets[3 * idx])

class TLE:
    """
    Two-line element set
    """
    def __init__(self, line1=None, line2=None, name=None):
        """
        Initialise TLE

        Parameters
        ----------
        line1, line2 : str, optional
            First and second lines of TLE
            Default = None (elements set to defaults)
        name : str, optional
            Name of object
            Default = None (name set to 'UNKNOWN')
        """
        if self._validate(line1, line2, name):
            self._unpack()
            self._ts = TS
            self._propagate = 0
        else:
            sys.exit()

    def _validate(self, line1, line2, name):
        """
        Validate TLE input

        Parameters
        ----------
        line1, line2 : str
            First and second lines of TLE
        name : str
            Name of object
        """
        self.line1 = line1
        self.line2 = line2
        for input in [self.line1, self.line2]:
            if input is not None:
                if not isinstance(input, str):
                    print('TLEError: Input must be str')
                    return 0
                if len(input) != 69:
                    print('TLEError: Lines must contain 69 characters')
                    return 0
        self.name = name
        if self.name is not None:
            if not isinstance(self.name, str):
                print('TLEError: Input must be str')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack TLE input
        """
        unpackElements(self)
        self._object = EarthSatellite(self.line1,
                                      self.line2,
                                      self.name)

    def _propagate_pos(self):
        """
        Predict geometric ICRS position(s) for stored time(s) and site

        Returns
        -------
        pos : skyfield.positionlib.Geometric
            Predicted geometric ICRS position(s)
        """
        return (self._object - self._site.topocentric).at(self._internal_time)
    
    def parse_propagation_info(self, time, site):
        """
        Parse time and site information for TLE propagation
        
        Parameters
        ----------
        time : datetime | array-like
            Input time | times for propagation
        site : str | st.site.Site
            Observation site name | Site object
        """
        if not parsePropagationInfo(self, time, site):
            sys.exit()

    def propagate_radec(self):
        """
        Predict sky coordinate(s) (RA, DEC) for stored time(s) and site

        Returns
        -------
        ra, dec : array-like | None
            Predicted RA [deg] and DEC [deg] | None if no time(s) or site stored
        """
        if self._propagate:
            ra, dec, _ = self._propagate_pos().radec()
            return ra._degrees, dec.degrees
        else:
            print('TLEError: No propagation time(s) or site stored')
            return None

    def propagate_ha(self):
        """
        Predict hour angle(s) for stored time(s) and site

        Returns
        -------
        ha : astropy.coordinates.Longitude | None
            Hour angle(s) [hours] | None if no time(s) or site stored
        """
        if self._propagate:
            ra, _ = self.propagate_radec()
            lst = Time(self._time, scale='utc', location=self._site.geodetic).sidereal_time('apparent')
            return (lst - ra).wrap_at(12 * u.hourangle)
        else:
            print('TLEError: No propagation time(s) or site stored')
            return None

    def propagate_altaz(self):
        """
        Predict sky coordinate(s) (ALT, AZ) for stored time(s) and site

        Returns
        -------
        coord : array-like | None
            Predicted ALT [deg] and AZ [deg] | None if no time(s) or site stored
        """
        if self._propagate:
            alt, az, _ = self._propagate_pos().altaz()
            return alt.degrees, az.degrees
        else:
            print('TLEError: No propagation time(s) or site stored')
            return None

    def modify_elements(self, checksum=False, norad_id=None, designator=None, epoch=None,
                        mmdot=None, mmdot2=None, drag=None, setnumber=None,
                        inclination=None, raan=None, eccentricity=None, argperigee=None,
                        meananomaly=None, mm=None, revnumber=None):
        """
        Modify elements of TLE

        Parameters
        ----------
        checksum : bool, optional
            Toggle to recalculate checksum
            Default = False (ignore checksum)

        See descriptions in st.tlelements for suitable inputs for kwargs, optional
        Defaults = None (no modification)
        """
        modifyElements(self,
                       checksum=checksum,
                       norad_id=norad_id,
                       designator=designator,
                       epoch=epoch,
                       mmdot=mmdot,
                       mmdot2=mmdot2,
                       drag=drag,
                       setnumber=setnumber,
                       inclination=inclination,
                       raan=raan,
                       eccentricity=eccentricity,
                       argperigee=argperigee,
                       meananomaly=meananomaly,
                       mm=mm,
                       revnumber=revnumber)
