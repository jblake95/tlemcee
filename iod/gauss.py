"""
Gauss method of preliminary orbit determination
"""

import numpy as np

from iod.utils import (
    Constants,
    positionVector,
    directionCosineVector,
    improveOrbit,
    orbitalElementsAlgorithm
)

def gaussAlgorithm(arc, mu=Constants.mu, f=Constants.f, r_body=Constants.R_e,
                   tolerance=1e-6, sig=5, improve=False, verbose=False):
    """
    Carry out Gauss method of preliminary orbit determination

    Parameters
    ----------
    arc : arc.Arc
        Orbital arc object containing ephemeris and site information
    mu : float, optional
        Gravitational parameter [km^3/s^2]
        Default = Constants.mu [Earth]
    f : float, optional
        Oblateness of central body
        Default = Constants.f [Earth]
    r_body : float, optional
        Radius of central body
        Default = Constants.R_e [Earth]
    tolerance : float, optional
        Tolerance for universal Kepler equation solver
        (used if improve=True)
        Default = 1e-6
    sig : int, optional
        Number of significant figures for improvement algorithm success
        (used if improve=True)
        Default = 5
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
        [1, i] for i orbital solutions, each containing:
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
    # obtain position and direction cosine vectors
    R_1 = positionVector(arc.site.geodetic.lat.rad,
                         arc.obs1.lst.rad,
                         arc.site.geodetic.height.value,
                         r_body,
                         f)
    R_2 = positionVector(arc.site.geodetic.lat.rad,
                         arc.obs2.lst.rad,
                         arc.site.geodetic.height.value,
                         r_body,
                         f)
    R_3 = positionVector(arc.site.geodetic.lat.rad,
                         arc.obs3.lst.rad,
                         arc.site.geodetic.height.value,
                         r_body,
                         f)

    rho_hat_1 = directionCosineVector(arc.obs1.coord.ra.rad,
                                      arc.obs1.coord.dec.rad)
    rho_hat_2 = directionCosineVector(arc.obs2.coord.ra.rad,
                                      arc.obs2.coord.dec.rad)
    rho_hat_3 = directionCosineVector(arc.obs3.coord.ra.rad,
                                      arc.obs3.coord.dec.rad)

    # step 1 - time intervals
    tau_1 = (arc.obs1.utc - arc.obs2.utc).total_seconds()
    tau_3 = (arc.obs3.utc - arc.obs2.utc).total_seconds()
    tau = tau_3 - tau_1

    output = {'tau_1': tau_1,
              'tau_3': tau_3,
              'tau': tau}

    # step 2 - rho_hat cross products
    p_1 = np.cross(rho_hat_2, rho_hat_3)
    p_2 = np.cross(rho_hat_1, rho_hat_3)
    p_3 = np.cross(rho_hat_1, rho_hat_2)

    # step 3 - rho_hat scalar triple product
    d_0 = np.dot(rho_hat_1, p_1)

    # step 4 - compute scalar quantities
    d_11 = np.dot(R_1, p_1)
    d_12 = np.dot(R_1, p_2)
    d_13 = np.dot(R_1, p_3)
    d_21 = np.dot(R_2, p_1)
    d_22 = np.dot(R_2, p_2)
    d_23 = np.dot(R_2, p_3)
    d_31 = np.dot(R_3, p_1)
    d_32 = np.dot(R_3, p_2)
    d_33 = np.dot(R_3, p_3)

    # step 5 - calculate scalar position coefficients
    A = (1 / d_0) * (-d_12 * (tau_3 / tau) + d_22 + d_32 * (tau_1 / tau))
    B = (1 / (6 * d_0)) * (d_12 * (tau_3 ** 2 - tau ** 2) * tau_3 / tau +
                           d_32 * (tau ** 2 - tau_1 ** 2) * tau_1 / tau)
    E = np.dot(R_2, rho_hat_2)

    # step 6 - squared scalar distance for obs2
    R2_2 = np.dot(R_2, R_2)

    # step 7 - coefficients of scalar distance polynomial for obs2
    a = -(A ** 2 + 2 * A * E + R2_2)
    b = -2 * mu * B * (A + E)
    c = -(mu ** 2) * (B ** 2)

    # step 8 - find root of scalar distance polynomial for obs2
    physical_roots = []
    for root in np.roots([1, 0, a, 0, 0, b, 0, 0, c]):
        if not np.iscomplex(root):
            if root >= 0:
                physical_roots.append(np.real(root))

    if len(physical_roots) == 0:
        print('warning: no physical roots - cannot solve.')
        return None
    else:
        for r, root in enumerate(physical_roots):  # possible geocentric radii r_2

            # step 9 - obtain slant ranges
            num_1 = (6 * (d_31 * (tau_1 / tau_3) + d_21 * (tau / tau_3)) * root ** 3 +
                     mu * d_31 * (tau ** 2 - tau_1 ** 2) * (tau_1 / tau_3))
            den_1 = 6 * root ** 3 + mu * (tau ** 2 - tau_3 ** 2)
            num_3 = (6 * (d_13 * (tau_3 / tau_1) - d_23 * (tau / tau_1)) * root ** 3 +
                     mu * d_13 * (tau ** 2 - tau_3 ** 2) * (tau_3 / tau_1))
            den_3 = 6 * root ** 3 + mu * (tau ** 2 - tau_1 ** 2)

            rho_1 = (1 / d_0) * (num_1 / den_1 - d_11)
            rho_2 = A + mu * B / root ** 3
            rho_3 = (1 / d_0) * (num_3 / den_3 - d_33)

            # step 10 - orbiting body geocentric position vectors
            r_1 = R_1 + rho_1 * rho_hat_1
            r_2 = R_2 + rho_2 * rho_hat_2
            r_3 = R_3 + rho_3 * rho_hat_3

            # step 11 - Lagrange coefficients
            f_1 = 1 - (1 / 2) * (mu / root ** 3) * tau_1 ** 2
            f_3 = 1 - (1 / 2) * (mu / root ** 3) * tau_3 ** 2
            g_1 = tau_1 - (1 / 6) * (mu / root ** 3) * tau_1 ** 3
            g_3 = tau_3 - (1 / 6) * (mu / root ** 3) * tau_3 ** 3

            # step 12 - velocity vector for obs2
            v_2 = (1 / (f_1 * g_3 - f_3 * g_1)) * (-f_3 * r_1 + f_1 * r_3)

            # step 13 - orbital elements
            orbit = orbitalElementsAlgorithm(r_2,
                                             v_2,
                                             mu,
                                             verbose=verbose)

            # Improve the state vector if requested
            if improve:
                params = {'r_vec': r_2, 'v_vec': v_2,
                          'tau_1': tau_1, 'tau_3': tau_3,
                          'R_1': R_1, 'R_2': R_2, 'R_3': R_3,
                          'rho_hat_1': rho_hat_1,
                          'rho_hat_2': rho_hat_2,
                          'rho_hat_3': rho_hat_3,
                          'rho_1': rho_1, 'rho_2': rho_2, 'rho_3': rho_3,
                          'f_1': f_1, 'f_3': f_3,
                          'g_1': g_1, 'g_3': g_3,
                          'd_0': d_0,
                          'd_11': d_11, 'd_12': d_12, 'd_13': d_13,
                          'd_21': d_21, 'd_22': d_22, 'd_23': d_23,
                          'd_31': d_31, 'd_32': d_32, 'd_33': d_33}
                improved_orbit = improveOrbit(params,
                                              mu,
                                              tolerance=tolerance,
                                              sig=sig,
                                              verbose=verbose)
                if improved_orbit is not None:
                    orbit = improved_orbit

            output.update({r: orbit})

    return output
