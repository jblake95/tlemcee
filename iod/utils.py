"""
General tasks for the Gauss method
"""

import numpy as np
from datetime import timedelta
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import Angle

from st.tle import ST

class Log:
    """
    Convenience class for log-keeping
    """
    def __init__(self, n_rows, col_names, col_dtype):
        """
        Initialise Log

        Parameters
        ----------
        n_rows : int
            Number of rows for log
        col_names : list
            List of column names, e.g. ['col0', 'col1']
        col_dtype : list
            List of column dtypes, e.g. ['U25', 'f8']
        """
        self.table = Table(np.zeros((n_rows, len(col_names))),
                           names=col_names,
                           dtype=col_dtype)

    def update(self, column, row, entry):
        """
        Update a cell within the log

        Parameters
        ----------
        column : str
            Name of column containing cell
        row : int
            Index of row containing cell
        entry : int, float, str
            Desired cell entry
        """
        self.table[column][row] = entry

    def add_row(self, entry):
        """
        Add a row to the log

        Parameters
        ----------
        entry : List
            List of cell entries to append to the log
        """
        self.table.add_row(entry)

    def save(self, out_path):
        """
        Save log to file

        Parameters
        ----------
        out_path : str
            Path to output file
        """
        self.table.write(out_path,
                         format='csv',
                         overwrite=True)

def simulateEphem(arc, norad_id, start, end, step, out_dir=None):
    """
    Simulate ephemeris information for an orbital arc

    Parameters
    ----------
    arc : arc.Arc
        Orbital arc containing site information
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
        Default=None

    Returns
    -------
    None
    """
    span = end - start
    assert span > timedelta(seconds=0), "start epoch must be before end epoch"

    tle = ST().get_past_tle(norad_id,
                            start + span / 2)

    times = []
    for i in range(int(span / step)):
        times.append(start + i * step)

    ephem = Log(len(times),
                ['UTC', 'RA', 'DEC', 'HA', 'ALT', 'AZ'],
                ['U30', 'f8', 'f8', 'f8', 'f8', 'f8'])

    tle.parse_propagation_info(times, arc.site)

    ra, dec = tle.propagate_radec()
    ha = tle.propagate_ha()
    alt, az = tle.propagate_altaz()

    for t, time in enumerate(times):
        time_entry = time.isoformat()
        if len(time_entry) > 30:
            time_entry = time_entry[:30]
        ephem.update('UTC', t, time_entry)
        ephem.update('RA', t, ra[t])
        ephem.update('DEC', t, dec[t])
        ephem.update('HA', t, ha[t].value)
        ephem.update('ALT', t, alt[t])
        ephem.update('AZ', t, az[t])

    if out_dir is not None:
        out_path = '{}/simulated_arc_{}_{}_{}.csv'.format(out_dir,
                                                          norad_id,
                                                          start.strftime('%Y%m%d'),
                                                          site.name)
        ephem.save(out_path)

    arc.ephem = ephem.table
    return None

class Constants:
    """
    Container for constants used in Gauss method algorithm
    """
    G = 6.6740831e-11  # gravitational constant [m^3/kg/s^2]
    M_e = 5.9722e+24   # Mass of Earth [kg]
    R_e = 6378000.     # Radius of Earth [m]
    f = 0.003353       # oblateness of Earth
    mu = 398600.       # gravitational parameter [km^3/s^2]

def positionVector(phi, lst, h, r_body, f):
    """
    Position vector of an observer at a given time

    Parameters
    ----------
    phi : float
        Geodetic latitude of observer's location - angle between
        equatorial and normal planes [rad]
    lst : float
        Local sidereal time for observation [rad]
    h : float
        Altitude of observer [m]
    r_body : float
        Radius of central body [m]
    f : float
        Oblateness of central body

    Returns
    -------
    r : array-like
        Position vector of observer for given time [km]
    """
    r_x = (r_body / np.sqrt(1 - (2 * f - f ** 2) * np.sin(phi) ** 2) + h) \
          * np.cos(phi) * np.cos(lst)
    r_y = (r_body / np.sqrt(1 - (2 * f - f ** 2) * np.sin(phi) ** 2) + h) \
          * np.cos(phi) * np.sin(lst)
    r_z = (r_body * (1 - f) ** 2 / np.sqrt(1 - (2 * f - f ** 2) * np.sin(phi) ** 2) + h) \
          * np.sin(phi)
    return np.array([r_x, r_y, r_z]) / 1000.

def directionCosineVector(alpha, delta):
    """
    Direction cosine vector for orbiting body

    Parameters
    ----------
    alpha, delta : float
        Topocentric right ascension & declination measured [rad]

    Returns
    -------
    rho_hat : array-like
        Direction cosine vector for orbiting body
    """
    rho_hat_x = np.cos(delta) * np.cos(alpha)
    rho_hat_y = np.cos(delta) * np.sin(alpha)
    rho_hat_z = np.sin(delta)
    return np.array([rho_hat_x, rho_hat_y, rho_hat_z])

def stumpffC(z):
    """
    Form of Stumpff function C(z) in terms of circular and hyperbolic
    trigonometric functions

    Parameters
    ----------
    z : float
        Input value, calculated from reciprocal semimajor axis (alpha) and
        universal anomaly (chi) as z = alpha * chi **2

    Returns
    -------
    c : float
        Output value
    """
    if z > 0:
        return (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        return (np.cosh(np.sqrt(-z)) - 1) / -z
    else:
        return 1 / 2

def stumpffS(z):
    """
    Form of Stumpff function S(z) in terms of circular and hyperbolic
    trigonometric functions

    Parameters
    ----------
    z : float
        Input value, calculated from reciprocal semimajor axis (alpha) and
        universal anomaly (chi) as z = alpha * chi **2

    Returns
    -------
    s : float
        Output value
    """
    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z) ** 3
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / np.sqrt(-z) ** 3
    else:
        return 1 / 6

def round_sig(x, sig=5):
    """
    Round a given number to requested number of significant figures

    Parameters
    ----------
    x : float
        Number to be rounded
    sig : int
        Number of significant figures

    Returns
    -------
    y : float
        Rounded number
    """
    if x == 0.:
        return 0.
    else:
        return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)

def solveUniversalKepler(delta_t, r_0, v_r0, alpha, mu, tolerance=1e-6):
    """
    Solve universal Kepler's equation for universal anomaly (chi)

    Parameters
    ----------
    delta_t : float
        Time interval
    r_0, v_r0 : float
        Orbital radius and radial velocity at time t = t_0
    alpha : float
        Reciprocal of semimajor axis
    mu : float
        Gravitational parameter [km^3/s^2]
    tolerance : float, optional
        Tolerance defining algorithmic success
        Default = 1e-6

    Returns
    -------
    chi : float | None
        Universal anomaly | None upon failure
    """
    # step 1 - initial estimate of chi
    chi = np.sqrt(mu) * abs(alpha) * delta_t

    count = 0
    ratio = np.inf
    while abs(ratio) > tolerance:
        if count > 10:
            print('warning: universal Kepler solver failed to achieve tolerance')
            return None

        # step 4 - if ratio (see below) exceeds tolerance, update chi
        if ratio is not np.inf:
            chi -= ratio

        # step 2 - calculate f(chi) & f'(chi)
        z = alpha * chi ** 2
        f_chi = ((r_0 * v_r0 / np.sqrt(mu)) * chi ** 2 * stumpffC(z) +
                 (1 - alpha * r_0) * chi ** 3 * stumpffS(z) + r_0 * chi -
                 np.sqrt(mu) * delta_t)
        f_chi_prime = ((r_0 * v_r0 / np.sqrt(mu)) * chi * (1 - z * stumpffS(z)) +
                       (1 - alpha * r_0) * chi ** 2 * stumpffC(z) + r_0)

        # step 3 - calculate ratio
        ratio = f_chi / f_chi_prime

        count += 1

    # step 5 - if ratio is less than tolerance, accept chi as solution
    return chi

def improveOrbit(p, mu, tolerance=1e-6, sig=5, verbose=False):
    """
    Iterative improvement of orbit determined with Gauss method

    Parameters
    ----------
    p : dict
        Dictionary of parameters carried forward from Gauss method
    mu : float
        Gravitational parameter [km^3/s^2]
    tolerance : float, optional
        Tolerance for universal Kepler equation solver
        Default = 1e-6
    sig : int, optional
        Number of significant figures for improvement algorithm success
        Default = 5
    verbose : bool, optional
        Toggle to display orbital information
        Default = False

    Returns
    -------
    orbit : dict | None
        Orbital information for orbiting body:
        a     - semimajor axis [km]
        r_p   - radius at perihelion [km]
        r_a   - radius at aphelion [km]
        e     - eccentricity
        i     - inclination [degrees]
        Omega - right ascension of the ascending node [degrees]
        omega - argument of perigee [degrees]
        theta - true anomaly [degrees]
        T     - period [hours],
        improved using 'exact' values of Lagrange coefficients | None upon failure
    """
    f_1_list = [p['f_1']]
    f_3_list = [p['f_3']]
    g_1_list = [p['g_1']]
    g_3_list = [p['g_3']]
    rho_1_list = [round_sig(p['rho_1'], sig=sig)]
    rho_2_list = [round_sig(p['rho_2'], sig=sig)]
    rho_3_list = [round_sig(p['rho_3'], sig=sig)]
    i = 0
    while True:
        # step 1 - distance and speed
        r = np.sqrt(np.dot(p['r_vec'], p['r_vec']))
        v = np.sqrt(np.dot(p['v_vec'], p['v_vec']))

        # step 2 - reciprocal of semimajor axis
        alpha = 2 / r - v ** 2 / mu

        # step 3 - radial component of v_vec
        v_r = np.dot(p['v_vec'], p['r_vec']) / r

        # step 4 - solve universal Kepler's equation for chi_i
        chi_1 = solveUniversalKepler(p['tau_1'],
                                     r,
                                     v_r,
                                     alpha,
                                     mu,
                                     tolerance=tolerance)
        chi_3 = solveUniversalKepler(p['tau_3'],
                                     r,
                                     v_r,
                                     alpha,
                                     mu,
                                     tolerance=tolerance)

        if (chi_1 is None) | (chi_3 is None):
            return None

        # step 5 - use chi_i to determine new Lagrange coefficients
        z_1 = alpha * chi_1 ** 2
        z_3 = alpha * chi_3 ** 2

        f_1 = 1 - (chi_1 ** 2 / r) * stumpffC(z_1)
        g_1 = p['tau_1'] - (1 / np.sqrt(mu)) * chi_1 ** 3 * stumpffS(z_1)
        f_3 = 1 - (chi_3 ** 2 / r) * stumpffC(z_3)
        g_3 = p['tau_3'] - (1 / np.sqrt(mu)) * chi_3 ** 3 * stumpffS(z_3)

        f_1_list.append(f_1)
        f_3_list.append(f_3)
        g_1_list.append(g_1)
        g_3_list.append(g_3)

        f_1 = (f_1_list[i] + f_1_list[i - 1]) / 2
        f_3 = (f_3_list[i] + f_3_list[i - 1]) / 2
        g_1 = (g_1_list[i] + g_1_list[i - 1]) / 2
        g_3 = (g_3_list[i] + g_3_list[i - 1]) / 2  # take average f_i & g_i

        # step 6 - determine c_i
        c_1 = g_3 / (f_1 * g_3 - f_3 * g_1)
        c_3 = -g_1 / (f_1 * g_3 - f_3 * g_1)

        # step 7 - update values of rho_i
        rho_1 = (1 / p['d_0']) * (-p['d_11'] + (1 / c_1) * p['d_21'] -
                                  (c_3 / c_1) * p['d_31'])
        rho_2 = (1 / p['d_0']) * (-c_1 * p['d_12'] + p['d_22'] -
                                  c_3 * p['d_32'])
        rho_3 = (1 / p['d_0']) * (-(c_1 / c_3) * p['d_13'] +
                                  (1 / c_3) * p['d_23'] - p['d_33'])

        rho_1_list.append(round_sig(rho_1, sig=sig))
        rho_2_list.append(round_sig(rho_2, sig=sig))
        rho_3_list.append(round_sig(rho_3, sig=sig))

        # step 8 - update geocentric position vectors r_i
        r_1 = p['R_1'] + rho_1 * p['rho_hat_1']
        r_2 = p['R_2'] + rho_2 * p['rho_hat_2']
        r_3 = p['R_3'] + rho_3 * p['rho_hat_3']

        # step 9 - update velocity vector v
        v_2 = (1 / (f_1 * g_3 - f_3 * g_1)) * (-f_3 * r_1 + f_1 * r_3)

        # step 10 - check if rho_i are the 'same' to within sig figs
        if (rho_1_list[i] == rho_1_list[i - 1] and
            rho_2_list[i] == rho_2_list[i - 1] and
            rho_3_list[i] == rho_3_list[i - 1]):
            break
        elif i > 100:
            print('warning: improvement algorithm failed to converge')
            return None
        else:
            p['r_vec'], p['v_vec'] = r_2, v_2  # update state vector and repeat

        i += 1

    return orbitalElementsAlgorithm(p['r_vec'],
                                    p['v_vec'],
                                    mu,
                                    verbose=verbose)

def orbitalElementsAlgorithm(r_vec, v_vec, mu, verbose=False):
    """
    Obtain orbital elements from state vector

    Parameters
    ----------
    r_vec : array-like
        Position vector for orbiting body [km]
    v_vec : array-like
        Velocity vector for orbiting body [km/s]
    mu : float
        Gravitational parameter [km^3/s^2]
    verbose : bool, optional
        Toggle to display orbital information
        Default = False

    Returns
    -------
    orbit : dict
        Orbital information for orbiting body:
        a     - semimajor axis [km]
        r_p   - radius at perihelion [km]
        r_a   - radius at aphelion [km]
        e     - eccentricity
        i     - inclination [degrees]
        Omega - right ascension of the ascending node [degrees]
        omega - argument of perigee [degrees]
        theta - true anomaly [degrees]
        T     - period [hours]
    """
    # step 1 - distance
    r = np.sqrt(np.dot(r_vec, r_vec))

    # step 2 - speed
    v = np.sqrt(np.dot(v_vec, v_vec))

    # step 3 - radial velocity
    v_r = np.dot(v_vec, r_vec) / r

    # step 4 - specific angular momentum
    h_vec = np.cross(r_vec, v_vec)

    # step 5 - magnitude of specific angular momentum
    h = np.sqrt(np.dot(h_vec, h_vec))

    # step 6 - inclination [rad]
    i = Angle(np.arccos(h_vec[2] / h),
              u.rad)

    # step 7 - node line vector
    k_hat = (0, 0, 1)
    N_vec = np.cross(k_hat, h_vec)

    # step 8 - magnitude of node line vector
    N = np.sqrt(np.dot(N_vec, N_vec))

    # step 9 - right ascension of the ascending node [rad]
    if N_vec[1] >= 0:
        Omega = Angle(np.arccos(N_vec[0] / N),
                      u.rad)
    else:
        Omega = Angle(2 * np.pi - np.arccos(N_vec[0] / N),
                      u.rad)

    # step 10 - eccentricity vector
    e_vec = (1 / mu) * ((v ** 2 - (mu / r)) * r_vec - r * v_r * v_vec)

    # step 11 - eccentricity
    e = np.sqrt(np.dot(e_vec, e_vec))

    # step 12 - argument of perigee
    if e_vec[2] >= 0:
        omega = Angle(np.arccos(np.dot(N_vec, e_vec) / (N * e)),
                      u.rad)
    else:
        omega = Angle(2 * np.pi - np.arccos(np.dot(N_vec, e_vec) / (N * e)),
                      u.rad)

    # step 13 - true anomaly
    if v_r >= 0:
        theta = Angle(np.arccos(np.dot(e_vec, r_vec) / (e * r)),
                      u.rad)
    else:
        theta = Angle(2 * np.pi - np.arccos(np.dot(e_vec, r_vec) / (e * r)),
                      u.rad)

    # extras - perigee & apogee radii, semimajor axis, period
    r_p = (h ** 2 / mu) * (1 / (1 + e * np.cos(0.)))
    r_a = (h ** 2 / mu) * (1 / (1 + e * np.cos(np.pi)))

    a = (1 / 2) * (r_p + r_a)
    T = (2 * np.pi / np.sqrt(mu)) * a ** (3 / 2)

    if verbose:
        print('--------------------\n'
              'Orbital information:\n'
              '--------------------\n'
              'a[km]      = {}\n'
              'e          = {}\n'
              'i[deg]     = {}\n'
              'Omega[deg] = {}\n'
              'omega[deg] = {}\n'
              'theta[deg] = {}\n'
              'T[hrs]     = {}\n'
              '--------------------'.format(a,
                                            e,
                                            i.deg,
                                            Omega.deg,
                                            omega.deg,
                                            theta.deg,
                                            T / 3600.))
    return {'a': a,
            'r_p': r_p,
            'r_a': r_a,
            'e': e,
            'i': i.deg,
            'Omega': Omega.deg,
            'omega': omega.deg,
            'theta': theta.deg,
            'T': T / 3600.}
