"""
Two-line element (TLE) manipulation
"""

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

        self._ts = TS
        self._object = EarthSatellite(self.line1,
                                      self.line2,
                                      self.name)

    def _specify_utc(self, t):
        """
        Correct timezone info for existing datetime

        Parameters
        ----------
        t : datetime.datetime
            Time to correct

        Returns
        -------
        t : datetime.datetime
            Corrected time
        """
        return t.replace(tzinfo=utc)

    def _pos_at(self, t, site):
        """
        Predict geometric ICRS position at time t

        Parameters
        ----------
        t : datetime.datetime
            Time [utc]
        site : site.Site
            Observation site

        Returns
        -------
        pos : skyfield.positionlib.Geometric
            Predicted geometric ICRS position at time t
        """
        return (self._object - site.topocentric).at(self._ts.utc(self._specify_utc(t)))

    def radec_at(self, t, site):
        """
        Predict sky coordinate(s) (RA, Dec) at time(s) t

        Parameters
        ----------
        t : datetime object
            Time(s) [utc]
        site : site.Site
            Observation site

        Returns
        -------
        coord : astropy.coordinates.SkyCoord
            Predicted sky coordinate(s) (RA, Dec) at time(s) t
        """
        ra, dec, _ = self._pos_at(t, site).radec()
        return SkyCoord(ra=ra.hours,
                        dec=dec.degrees,
                        unit=(u.hourangle, u.deg),
                        frame='icrs')

    def hourangle_at(self, t, site):
        """
        Predict hour angle(s) at time(s) t

        Parameters
        ----------
        t : datetime object
            Time(s) [utc]
        site : site.Site
            Observation site

        Returns
        -------
        ha : Longitude object
            Hour angle(s) of object at time(s) t
        """
        ra = self.radec_at(t, site).ra
        lst = Time(t, scale='utc', location=site.geodetic).sidereal_time('apparent')
        return (lst - ra).wrap_at(12 * u.hourangle)

    def altaz_at(self, t, site):
        """
        Predict sky coordinate(s) (AltAz frame) at time(s) t

        Parameters
        ----------
        t : datetime object
            Time(s) [utc]
        site : site.Site
            Observation site

        Returns
        -------
        coord : astropy.coordinates.SkyCoord
            Predicted sky coordinate(s) (AltAz frame) at time(s) t
        """
        alt, az, _ = self._pos_at(t, site).altaz()
        return SkyCoord(alt=alt.degrees,
                        az=az.degrees,
                        unit=(u.deg, u.deg),
                        frame='altaz')

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
