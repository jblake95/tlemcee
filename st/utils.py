"""
General tasks for TLE manipulation
"""

def fractional_yearday(dt):
    """
    Obtain fractional year day from datetime object

    Parameters
    ----------
    dt : datetime.datetime
        datetime containing epoch information

    Returns
    -------
    year_day : float
        Fractional year day
    """
    frac = dt.hour / 24. + \
           dt.minute / (24. * 60.) + \
           dt.second / (24. * 60. * 60.) + \
           dt.microsecond / (24. * 60. * 60. * 1e6)

    return round(dt.timetuple().tm_yday + frac, 8)

def n_digits(integer):
    """
    Determine number of digits in a given integer

    Parameters
    ----------
    integer : int
        Input integer

    Returns
    -------
    n_digits : int
        Number of digits in input integer
    """
    return len(str(abs(int(integer))))

def tle_standard_form(number):
    """
    Convert number to scientific notation string consistent with TLE format

    Parameters
    ----------
    number : float
        Number to convert to TLE-friendly scientific notation

    Returns
    -------
    sci_not : str
        Correctly formatted string giving TLE-friendly scientific
        notation for input number
    """
    sci_not = '{:.4e}'.format(number)
    e_pos = sci_not.find('e')

    mantissa = float(sci_not[:e_pos])
    if mantissa >= 0:
        mantissa_sign = '+'
    else:
        mantissa_sign = '-'
    mantissa_str = '{:.5f}'.format(mantissa / 10.).lstrip('+-0')
    mantissa_str = mantissa_str.lstrip('.')

    if mantissa != 0.:
        exponent = int(sci_not[e_pos + 1:]) + 1
    else:
        exponent = 0
    if exponent >= 0:
        exponent_sign = '+'
    else:
        exponent_sign = '-'
    exponent_str = str(exponent).lstrip('+-')

    return '{}{}{}{}'.format(mantissa_sign,
                             mantissa_str,
                             exponent_sign,
                             exponent_str)

def read_tle_standard_form(sci_not):
    """
    Convert TLE scientific notation string to number

    Parameters
    ----------
    sci_not : str
        TLE scientific notation string
    """
    sci_not = '.{}'.format(sci_not.lstrip('+- '))

    return float('{}e{}'.format(sci_not[:-2],
                                sci_not[-2:]))

