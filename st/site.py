"""
Class for observatory information
"""

from skyfield.api import Topos
from astropy.coordinates import EarthLocation

# add new observation sites here
# format {site_name: [longitude (deg), latitude (deg), elevation (m)]}
SITES = {'LaPalma': [-17.877594, 28.761938, 2348]}

class Site:
    """
    Observation site
    """
    def __init__(self, site_name):
        """
        Initialise Site

        Parameters
        ----------
        site_name : str
            Identifier for desired observation site, e.g. 'LaPalma'
        """
        if self._validate(site_name):
            self._set_location()

    def _validate(self, site_name):
        """
        Validate input

        Parameters
        ----------
        site_name : str
            Identifier for desired observation site, e.g. 'LaPalma'
        """
        self.name = site_name
        if not isinstance(self.name, str):
            print('SiteError: name must be string')
            return 0
        if not {self.name} <= set(SITES):
            print('SiteError: name not recognised')
            return 0
        return 1

    def _set_location(self):
        """
        Set site location
        """
        self.longitude = SITES[self.name][0]
        self.latitude = SITES[self.name][1]
        self.elevation = SITES[self.name][2]

        self.geodetic = EarthLocation(lon=self.longitude,
                                      lat=self.latitude,
                                      height=self.elevation)

        self.topocentric = Topos(longitude_degrees=self.longitude,
                                 latitude_degrees=self.latitude,
                                 elevation_m=self.elevation)
