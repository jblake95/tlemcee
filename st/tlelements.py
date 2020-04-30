"""
Classes for TLE components
"""

import sys
from datetime import (
    datetime,
    timedelta
)
from astropy import units as u
from astropy.coordinates import Angle
from skyfield.sgp4lib import EarthSatellite

from st.utils import (
    fractional_yearday,
    n_digits,
    tle_standard_form,
    read_tle_standard_form
)

# line 1 components
class Norad_ID:
    """
    NORAD identification number
    """
    def __init__(self, norad_id=0):
        """
        Initialise Norad_ID

        Parameters
        ----------
        norad_id : str | int, optional
            TLE input (line 1) | value input
            Default = 0
        """
        if self._validate(norad_id):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Norad_ID
        """
        return 'Norad_ID <{}>'.format(self.entry)

    def _validate(self, norad_id):
        """
        Validate input

        Parameters
        ----------
        norad_id : str | int
            TLE input (line 1) | value input
        """
        self._input = norad_id
        if not isinstance(self._input, (str, int)):
            print('Norad_IDError: Input must be str | int')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('EpochError: TLE line 1 must contain 69 characters')
                return 0
        else:
            if self._input < 0:
                print('Norad_IDError: Input must be positive')
                return 0
            if n_digits(self._input) > 5:
                print('Norad_IDError: Input must not exceed 5 digits')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load NORAD identification number from TLE input
        """
        self.value = int(self._input[2:7])

    def _from_value(self):
        """
        Load NORAD identification number from value input
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = '{:5.0f}'.format(self.value)

class Designator:
    """
    International designator
    """
    def __init__(self, designator=None):
        """
        Initialise Designator

        Parameters
        ----------
        designator : str | dict, optional
            TLE input (line 1) | dict input with format:
            {
             'yr': Launch year, e.g. 2000 (int),
             'no': Launch number, e.g. 1 (int),
             'id': Launch ID, e.g. 'A' (str)
             }
            Default = None - empty (no designator)
        """
        if self._validate(designator):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Launch
        """
        return 'International designator <{}>'.format(self.entry)

    def _validate(self, designator):
        """
        Validate input

        Parameters
        ----------
        designator : str | dict
            TLE input (line 1) | dict input [key (value): 'yr' (int), 'no' (int), 'id' (str)]
        """
        self._input = designator
        if (self._input is not None) & (not isinstance(self._input, (str, dict))):
            print('DesignatorError: Input must be str | dict')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('DesignatorError: TLE line 1 must contain 69 characters')
                return 0

        if isinstance(self._input, dict):
            if not {'yr', 'no', 'id'} <= set(self._input):
                print('DesignatorError: dict input must contain keys: "yr", "no", "id"')
                return 0
            if not isinstance(self._input['yr'], int):
                print('DesignatorError: Input launch year must be int')
                return 0
            if isinstance(self._input['no'], int):
                if n_digits(self._input['no']) > 3:
                    print('DesignatorError: Input launch number exceeds limit')
                    return 0
            else:
                print('DesignatorError: Input launch number must be int')
                return 0
            if isinstance(self._input['id'], str):
                if len(self._input['id']) > 3:
                    print('DesignatorError: Input launch ID exceeds limit')
                    return 0
                if not self._input.isalpha():
                    print('DesignatorError: Invalid launch ID: ABC not 123')
            else:
                print('DesignatorError: Input launch ID must be str')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if self._input is not None:
            if isinstance(self._input, str):
                self._from_tle()
            else:
                self._from_dict()

        self._format_tle()

    def _from_tle(self):
        """
        Load Designator information from TLE input
        """
        year = self._input[9:11]
        if int(year) < 57:
            year = '20{}'.format(year)
        else:
            year = '19{}'.format(year)

        self.year = int(year)
        self.number = int(self._input[11:14])
        self.id = self._input[14:17].strip()

    def _from_dict(self):
        """
        Load Designator information from dict of components
        """
        self.year = self._input['yr']
        self.number = self._input['no']
        self.id = self._input['id']

    def _format_tle(self):
        """
        TLE-friendly format
        """
        if self._input is not None:
            self.entry = '{:02.0f}{:03.0f}{:3}'.format(self.year % 100,
                                                       self.number,
                                                       self.id)
        else:
            self.entry = '{:8}'.format('')

class Epoch:
    """
    TLE reference epoch
    """
    def __init__(self, epoch=None):
        """
        Initialise Epoch

        Parameters
        ----------
        epoch : str | datetime.datetime, optional
            TLE input (line 1) | datetime input
            Default = None - 01/01/2000 (warning: TLEs quickly become out-of-date!)
        """
        # validate input and load Epoch information
        if self._validate(epoch):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Epoch
        """
        return 'Reference epoch <{}>'.format(self.entry)

    def _validate(self, epoch):
        """
        Validate input

        Parameters
        ----------
        epoch : str | datetime.datetime
            TLE input (line 1) | datetime input
        """
        self._input = epoch
        if (self._input is not None) & (not isinstance(self._input, (str, datetime))):
            print('EpochError: Input must be str | datetime')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('EpochError: TLE line 1 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if self._input is not None:
            if isinstance(self._input, str):
                self._from_tle()
            else:
                self._from_datetime()
        else:
            self._from_default()

        self.year = self.date.year
        self.yearday = fractional_yearday(self.date)

        self._format_tle()

    def _from_tle(self):
        """
        Load Epoch information from TLE input
        """
        yearday = float(self._input[20:32])
        year = self._input[18:20]

        if int(year) < 57:
            year = '20{}'.format(year)
        else:
            year = '19{}'.format(year)

        self.date = datetime.strptime(year, '%Y')
        self.date += timedelta(days=yearday - 1.)

    def _from_datetime(self):
        """
        Load Epoch information from datetime input
        """
        self.date = self._input

    def _from_default(self):
        """
        Load Epoch information from default
        """
        self.date = datetime.strptime('2000', '%Y')

    def _format_tle(self):
        """
        TLE-friendly format
        """
        year_tle = '{:02.0f}'.format(self.year % 100)
        yearday_tle = '{:012.8f}'.format(self.yearday)

        self.entry = '{}{}'.format(year_tle,
                                   yearday_tle)

class Mmdot:
    """
    First time derivative of mean motion
    """
    def __init__(self, mmdot=0.):
        """
        Initialise Mmdot

        Parameters
        ----------
        mmdot : str | float, optional
            TLE input (line 1) | value input
            Default = 0.
        """
        if self._validate(mmdot):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Mmdot
        """
        return 'Mmdot <{}>'.format(self.entry)

    def _validate(self, mmdot):
        """
        Validate input

        Parameters
        ----------
        mmdot : str | float
            First time derivative of mean motion
        """
        self._input = mmdot
        if not isinstance(self._input, (str, float)):
            print('MmdotError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('Mmdot2Error: TLE line 1 must contain 69 characters')
                return 0
        else:
            if abs(mmdot) > 1:
                print('MmdotError: Input value exceeds limit')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load Mmdot information from TLE input
        """
        self.value = float(self._input[33:43])

    def _from_value(self):
        """
        Load Mmdot information from value input
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        if self.value >= 0:
            self.entry = '+' + '{:9.8f}'.format(self.value).lstrip('0')
        else:
            self.entry = '-' + '{:9.8f}'.format(self.value).lstrip('-0')

class Mmdot2:
    """
    Second time derivative of mean motion
    """
    def __init__(self, mmdot2=0.):
        """
        Initialise Mmdot2

        Parameters
        ----------
        mmdot2 : str | float, optional
            TLE input (line 1) | value input
            Default = 0.
        """
        if self._validate(mmdot2):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Mmdot2
        """
        return 'Mmdot2 <{}>'.format(self.entry)

    def _validate(self, mmdot2):
        """
        Validate input

        Parameters
        ----------
        mmdot2 : str | float
            TLE input (line 1) | value input
        """
        self._input = mmdot2
        if not isinstance(self._input, (str, float)):
            print('Mmdot2Error: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('Mmdot2Error: TLE line 1 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load Mmdot2 information from TLE input
        """
        self.value = read_tle_standard_form(self._input[44:52])

    def _from_value(self):
        """
        Load Mmdot2 information from value input
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = tle_standard_form(self.value)

class Drag:
    """
    B-star drag term
    """
    def __init__(self, drag=0.):
        """
        Initialise Drag

        Parameters
        ----------
        drag : str | float, optional
            TLE input (line 1) | value input
            Default = 0.
        """
        if self._validate(drag):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Drag
        """
        return 'Drag <{}>'.format(self.entry)

    def _validate(self, drag):
        """
        Validate input

        Parameters
        ----------
        drag : str | float
            TLE input (line 1) | value input
        """
        self._input = drag
        if not isinstance(self._input, (str, float)):
            print('DragError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('DragError: TLE line 1 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load Drag information from TLE input
        """
        self.value = read_tle_standard_form(self._input[53:61])

    def _from_value(self):
        """
        Load Drag information from value input
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = tle_standard_form(self.value)

class SetNumber:
    """
    Element set number
    """
    def __init__(self, set_no=0):
        """
        Initialise SetNumber

        Parameters
        ----------
        set_no : str | int, optional
            TLE input (line 1) | value input
            Default = 0
        """
        if self._validate(set_no):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for SetNumber
        """
        return 'SetNumber <{}>'.format(self.entry)

    def _validate(self, set_no):
        """
        Validate input

        Parameters
        ----------
        set_no : str | int
            TLE input (line 1) | value input
        """
        self._input = set_no
        if not isinstance(self._input, (str, int)):
            print('SetNumberError: Input must be str | int')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('SetNumberError: TLE line 1 must contain 69 characters')
                return 0
        else:
            if n_digits(self._input) > 4:
                print('SetNumberError: Input element set number exceeds limit')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load SetNumber information from TLE input
        """
        self.value = int(self._input[64:68])

    def _from_value(self):
        """
        Load SetNumber information from value input
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = '{:4}'.format(self.value)

class CheckSum:
    """
    TLE line checksum
    """
    def __init__(self, line):
        """
        Initialise CheckSum

        Parameters
        ----------
        line : str
            TLE input (line 1 or 2)
        """
        if self._validate(line):
            self._calculate()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for CheckSum
        """
        return 'CheckSum <{}>'.format(self.entry)

    def _validate(self, line):
        """
        Validate input

        Parameters
        ----------
        line : str
            TLE input (line 1 or 2)
        """
        self._input = line
        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('CheckSumError: TLE line must contain 69 characters')
                return 0
        else:
            print('CheckSumError: Input must be str')
            return 0
        return 1

    def _calculate(self):
        """
        Calculate checksum for line
        """
        self.value = 0
        for character in self._input[:-1]:
            if character.isnumeric():
                self.value += int(character)
            elif character == '-':
                self.value += 1
            else:
                continue
        self.value %= 10
        self._format_tle()

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = '{:1}'.format(self.value)

class Inclination:
    """
    Orbital inclination
    """
    def __init__(self, inclination=0.):
        """
        Initialise Inclination

        Parameters
        ----------
        inclination : str | float, optional
            TLE input (line 2) | value input [deg]
            Default = 0.
        """
        if self._validate(inclination):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Inclination
        """
        return 'Inclination <{}>'.format(self.entry)

    def _validate(self, inclination):
        """
        Validate input

        Parameters
        ----------
        inclination : str | float
            TLE input (line 2) | value input [deg]
        """
        self._input = inclination
        if not isinstance(self._input, (str, float)):
            print('InclinationError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('InclinationError: TLE line must contain 69 characters')
                return 0
        else:
            if not Angle(self._input, u.deg).is_within_bounds(0. * u.deg, 180. * u.deg):
                print('InclinationError: Input is outside interval [0, 180]')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load Inclination information from TLE input
        """
        self.value = Angle(float(self._input[8:16]), u.deg)

    def _from_value(self):
        """
        Load Inclination information from value input [deg]
        """
        self.value = Angle(self._input, u.deg)

    def _format_tle(self):
        """
        TLE-firendly format
        """
        self.entry = '{:8.4f}'.format(self.value.deg)

class RAAN:
    """
    Right ascension of ascending node
    """
    def __init__(self, raan=0.):
        """
        Initialise RAAN

        Parameters
        ----------
        raan : str | float, optional
            TLE input (line 2) | value input [deg]
            Default = 0.
        """
        if self._validate(raan):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for RAAN
        """
        return 'RAAN <{}>'.format(self.entry)

    def _validate(self, raan):
        """
        Validate input

        Parameters
        ----------
        raan : str | float
            TLE input (line 2) | value input [deg]
        """
        self._input = raan
        if not isinstance(self._input, (str, float)):
            print('RAANError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('RAANError: TLE line 2 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load RAAN information from TLE input
        """
        self.value = Angle(float(self._input[17:25]), u.deg).wrap_at(360. * u.deg)

    def _from_value(self):
        """
        Load RAAN information from value input [deg]
        """
        self.value = Angle(self._input, u.deg).wrap_at(360. * u.deg)

    def _format_tle(self):
        """
        TLE-firendly format
        """
        self.entry = '{:8.4f}'.format(self.value.deg)

class Eccentricity:
    """
    Orbital eccentricity
    """
    def __init__(self, eccentricity=0.):
        """
        Inititialise Eccentricity

        Parameters
        ----------
        eccentricity : str | float, optional
            TLE input (line 2) | value input
            Default = 0.
        """
        if self._validate(eccentricity):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Eccentricity
        """
        return 'Eccentricity <{}>'.format(self.entry)

    def _validate(self, eccentricity):
        """
        Validate input

        Parameters
        ----------
        eccentricity : str | float
            TLE input (line 2) | value input
        """
        self._input = eccentricity
        if not isinstance(self._input, (str, float)):
            print('EccentricityError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('EccentricityError: TLE line 2 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load Eccentricity information from TLE input
        """
        self.value = float('.{}'.format(self._input[26:33]))

    def _from_value(self):
        """
        Load Eccentricity information from value input
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = '{:.7f}'.format(self.value).split('.')[1]

class ArgPerigee:
    """
    Argument of perigee
    """
    def __init__(self, argp=0.):
        """
        Initialise ArgPerigee

        Parameters
        ----------
        argp : str | float, optional
            TLE input (line 2) | value input [deg]
            Default = 0.
        """
        if self._validate(argp):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for ArgPerigee
        """
        return 'ArgPerigee <{}>'.format(self.entry)

    def _validate(self, argp):
        """
        Validate input

        Parameters
        ----------
        argp : str | float
            TLE input (line 2) | value input [deg]
        """
        self._input = argp
        if not isinstance(self._input, (str, float)):
            print('ArgPerigeeError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('ArgPerigeeError: TLE line 2 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load ArgPerigee information from TLE input
        """
        self.value = Angle(float(self._input[34:42]), u.deg).wrap_at(360. * u.deg)

    def _from_value(self):
        """
        Load ArgPerigee information from value input [deg]
        """
        self.value = Angle(self._input, u.deg).wrap_at(360. * u.deg)

    def _format_tle(self):
        """
        TLE-firendly format
        """
        self.entry = '{:8.4f}'.format(self.value.deg)

class MeanAnomaly:
    """
    Mean anomaly
    """
    def __init__(self, anomaly=0.):
        """
        Initialise MeanAnomaly

        Parameters
        ----------
        anomaly : str | float, optional
            TLE input (line 2) | value input
            Default = 0.
        """
        if self._validate(anomaly):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for MeanAnomaly
        """
        return 'MeanAnomaly <{}>'.format(self.entry)

    def _validate(self, anomaly):
        """
        Validate input

        Parameters
        ----------
        anomaly : str | float
            TLE input (line 2) | value input [deg]
        """
        self._input = anomaly
        if not isinstance(self._input, (str, float)):
            print('MeanAnomalyError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('MeanAnomalyError: TLE line 2 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load MeanAnomaly information from TLE input
        """
        self.value = Angle(float(self._input[43:51]), u.deg).wrap_at(360. * u.deg)

    def _from_value(self):
        """
        Load MeanAnomaly information from value input [deg]
        """
        self.value = Angle(self._input, u.deg).wrap_at(360. * u.deg)

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = '{:8.4f}'.format(self.value.deg)

class Mm:
    """
    Mean motion
    """
    def __init__(self, mm=0.):
        """
        Initialise Mm

        Parameters
        ----------
        mm : str | float, optional
            TLE input (line 2) | value input [revs per day]
            Default = 0.
        """
        if self._validate(mm):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for Mm
        """
        return 'Mm <{}>'.format(self.entry)

    def _validate(self, mm):
        """
        Validate input

        Parameters
        ----------
        mm : str | float
            TLE input (line 2) | value input
        """
        self._input = mm
        if not isinstance(self._input, (str, float)):
            print('MmError: Input must be str | float')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('MmError: TLE line 2 must contain 69 characters')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load Mm information from TLE input
        """
        self.value = float(self._input[52:63])

    def _from_value(self):
        """
        Load Mm information from value input [revs per day]
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = '{:11.8f}'.format(self.value)

class RevNumber:
    """
    Revolution number
    """
    def __init__(self, rev_no=0):
        """
        Initialise RevNumber

        Parameters
        ----------
        rev_no : str | int, optional
            TLE input (line 2) | value input
            Default = 0
        """
        if self._validate(rev_no):
            self._unpack()
        else:
            sys.exit()

    def __str__(self):
        """
        Printable string for RevNumber
        """
        return 'SetNumber <{}>'.format(self.entry)

    def _validate(self, rev_no):
        """
        Validate input

        Parameters
        ----------
        rev_no : str | int
            TLE input (line 2) | value input
        """
        self._input = rev_no
        if not isinstance(self._input, (str, int)):
            print('RevNumberError: Input must be str | int')
            return 0

        if isinstance(self._input, str):
            if len(self._input) != 69:
                print('RevNumberError: TLE line 2 must contain 69 characters')
                return 0
        else:
            if n_digits(self._input) > 5:
                print('RevNumberError: Input revolution number exceeds limit')
                return 0
        return 1

    def _unpack(self):
        """
        Unpack input
        """
        if isinstance(self._input, str):
            self._from_tle()
        else:
            self._from_value()

        self._format_tle()

    def _from_tle(self):
        """
        Load RevNumber information from TLE input
        """
        self.value = int(self._input[63:68])

    def _from_value(self):
        """
        Load RevNumber information from value input
        """
        self.value = self._input

    def _format_tle(self):
        """
        TLE-friendly format
        """
        self.entry = '{:5}'.format(self.value)

class EmptyTLE:
    """
    Container for empty TLE lines
    """
    line1 = '1 {}U {} {} {} {} {} 0 {}0'
    line2 = '2 {} {} {} {} {} {} {}{}0'

def unpackElements(tle):
    """
    Unpack TLE input

    Parameters
    ----------
    tle : st.tle.TLE
        TLE object to unpack

    Returns
    -------
    None
    """
    tle.checksums = []
    # line 1 elements
    if tle.line1 is not None:
        tle.norad_id = Norad_ID(tle.line1)
        tle.designator = Designator(tle.line1)
        tle.epoch = Epoch(tle.line1)
        tle.mmdot = Mmdot(tle.line1)
        tle.mmdot2 = Mmdot2(tle.line1)
        tle.drag = Drag(tle.line1)
        tle.setnumber = SetNumber(tle.line1)
        tle.checksums.append(CheckSum(tle.line1))
    else:
        tle.norad_id = Norad_ID()
        tle.designator = Designator()
        tle.epoch = Epoch()
        tle.mmdot = Mmdot()
        tle.mmdot2 = Mmdot2()
        tle.drag = Drag()
        tle.setnumber = SetNumber()

        tle.line1 = EmptyTLE.line1.format(tle.norad_id.entry,
                                          tle.designator.entry,
                                          tle.epoch.entry,
                                          tle.mmdot.entry,
                                          tle.mmdot2.entry,
                                          tle.drag.entry,
                                          tle.setnumber.entry)
        # compute checksum
        tle.checksums.append(CheckSum(tle.line1))
        tle.line1 = '{}{}'.format(tle.line1[:-1],
                                  tle.checksums[0].entry)

    # line 2 elements
    if tle.line2 is not None:
        tle.inclination = Inclination(tle.line2)
        tle.raan = RAAN(tle.line2)
        tle.eccentricity = Eccentricity(tle.line2)
        tle.argperigee = ArgPerigee(tle.line2)
        tle.meananomaly = MeanAnomaly(tle.line2)
        tle.mm = Mm(tle.line2)
        tle.revnumber = RevNumber(tle.line2)
        tle.checksums.append(CheckSum(tle.line2))
    else:
        tle.inclination = Inclination()
        tle.raan = RAAN()
        tle.eccentricity = Eccentricity()
        tle.argperigee = ArgPerigee()
        tle.meananomaly = MeanAnomaly()
        tle.mm = Mm()
        tle.revnumber = RevNumber()

        tle.line2 = EmptyTLE.line2.format(tle.norad_id.entry,
                                          tle.inclination.entry,
                                          tle.raan.entry,
                                          tle.eccentricity.entry,
                                          tle.argperigee.entry,
                                          tle.meananomaly.entry,
                                          tle.mm.entry,
                                          tle.revnumber.entry)
        # compute checksum
        tle.checksums.append(CheckSum(tle.line2))
        tle.line2 = '{}{}'.format(tle.line2[:-1],
                                  tle.checksums[1].entry)

    # name
    if tle.name is not None:
        # remove line number if from 3le
        if tle.name[0] == '0':
            tle.name = tle.name[2:]
    else:
        tle.name = 'UNKNOWN'

    tle._object = EarthSatellite(tle.line1,
                                 tle.line2,
                                 tle.name)
    return None

def calculateCheckSums(tle):
    """
    Calculate checksums for TLE

    Parameters
    ----------
    tle : st.tle.TLE
        TLE object in need of checksum calculation

    Returns
    -------
    None
    """
    tle.checksums[0] = CheckSum(tle.line1)
    tle.line1 = '{}{}'.format(tle.line1[:-1],
                              tle.checksums[0].entry)
    tle.checksums[1] = CheckSum(tle.line2)
    tle.line2 = '{}{}'.format(tle.line2[:-1],
                              tle.checksums[1].entry)

def modifyElements(tle, checksum, norad_id, designator, epoch, mmdot, mmdot2, drag, setnumber,
                   inclination, raan, eccentricity, argperigee, meananomaly, mm, revnumber):
    """
    Modify elements of a TLE object

    Parameters
    ----------
    tle : st.tle.TLE
        TLE object in need of modification
    checksum : bool
        Toggle to recalculate checksum

    See descriptions above for suitable inputs for args
    (set to None for no modification)

    Returns
    -------
    None
    """
    if norad_id is not None:
        tle.norad_id = Norad_ID(norad_id)
        tle.line1 = '{}{}{}'.format(tle.line1[:2],
                                    tle.norad_id.entry,
                                    tle.line1[7:])
        tle.line2 = '{}{}{}'.format(tle.line2[:2],
                                    tle.norad_id.entry,
                                    tle.line2[7:])

    if designator is not None:
        tle.designator = Designator(designator)
        tle.line1 = '{}{}{}'.format(tle.line1[:9],
                                    tle.designator.entry,
                                    tle.line1[17:])

    if epoch is not None:
        tle.epoch = Epoch(epoch)
        tle.line1 = '{}{}{}'.format(tle.line1[:18],
                                    tle.epoch.entry,
                                    tle.line1[32:])

    if mmdot is not None:
        tle.mmdot = Mmdot(mmdot)
        tle.line1 = '{}{}{}'.format(tle.line1[:33],
                                    tle.mmdot.entry,
                                    tle.line1[43:])

    if mmdot2 is not None:
        tle.mmdot2 = Mmdot2(mmdot2)
        tle.line1 = '{}{}{}'.format(tle.line1[:44],
                                    tle.mmdot2.entry,
                                    tle.line1[52:])

    if drag is not None:
        tle.drag = Drag(drag)
        tle.line1 = '{}{}{}'.format(tle.line1[:53],
                                    tle.drag.entry,
                                    tle.line1[61:])

    if setnumber is not None:
        tle.setnumber = SetNumber(setnumber)
        tle.line1 = '{}{}{}'.format(tle.line1[:64],
                                    tle.setnumber.entry,
                                    tle.line1[68:])

    if inclination is not None:
        tle.inclination = Inclination(inclination)
        tle.line2 = '{}{}{}'.format(tle.line2[:8],
                                    tle.inclination.entry,
                                    tle.line2[16:])

    if raan is not None:
        tle.raan = RAAN(raan)
        tle.line2 = '{}{}{}'.format(tle.line2[:17],
                                    tle.raan.entry,
                                    tle.line2[25:])

    if eccentricity is not None:
        tle.eccentricity = Eccentricity(eccentricity)
        tle.line2 = '{}{}{}'.format(tle.line2[:26],
                                    tle.eccentricity.entry,
                                    tle.line2[33:])

    if argperigee is not None:
        tle.argperigee = ArgPerigee(argperigee)
        tle.line2 = '{}{}{}'.format(tle.line2[:34],
                                    tle.argperigee.entry,
                                    tle.line2[42:])

    if meananomaly is not None:
        tle.meananomaly = MeanAnomaly(meananomaly)
        tle.line2 = '{}{}{}'.format(tle.line2[:43],
                                    tle.meananomaly.entry,
                                    tle.line2[51:])

    if mm is not None:
        tle.mm = Mm(mm)
        tle.line2 = '{}{}{}'.format(tle.line2[:52],
                                    tle.mm.entry,
                                    tle.line2[63:])

    if revnumber is not None:
        tle.revnumber = RevNumber(revnumber)
        tle.line2 = '{}{}{}'.format(tle.line2[:63],
                                    tle.revnumber.entry,
                                    tle.line2[68:])

    if checksum:
        calculateCheckSums(tle)

    tle._object = EarthSatellite(tle.line1,
                                 tle.line2,
                                 tle.name)
    return None
