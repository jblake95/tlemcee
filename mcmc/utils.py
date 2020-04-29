"""
General tasks for tlemcee MCMC algorithm
"""
from datetime import datetime
import emcee
import numpy as np
import time

from mcmc.config import Config

def model(theta, config):
    """
    Model: modify elements and compute new arc

    Parameters
    ----------
    theta : array-like
        Current set of parameters from MCMC
    config : mcmc.Config
        Config containing prior and data information

    Returns
    -------
    ra, dec : array-like
        Updated arc, computed with modified TLE
    """
    t0 = time.time()
    params = config.priors['name']
    if 'mmdot' in params:
        mmdot = float(theta[params.index('mmdot')])
    else:
        mmdot = None
    if 'mmdot2' in params:
        mmdot2 = float(theta[params.index('mmdot2')])
    else:
        mmdot2 = None
    if 'drag' in params:
        drag = float(theta[params.index('drag')])
    else:
        drag = None
    if 'inclination' in params:
        inclination = float(theta[params.index('inclination')])
    else:
        inclination = None
    if 'raan' in params:
        raan = float(theta[params.index('raan')])
    else:
        raan = None
    if 'eccentricity' in params:
        eccentricity = float(theta[params.index('eccentricity')])
    else:
        eccentricity = None
    if 'argperigee' in params:
        argperigee = float(theta[params.index('argperigee')])
    else:
        argperigee = None
    if 'meananomaly' in params:
        meananomaly = float(theta[params.index('meananomaly')])
    else:
        meananomaly = None
    if 'mm' in params:
        mm = float(theta[params.index('mm')])
    else:
        mm = None
    t1 = time.time()
    print(t1 - t0)
    t0 = time.time()
    config.tle.modify_elements(mmdot=mmdot,
                               mmdot2=mmdot2,
                               drag=drag,
                               inclination=inclination,
                               raan=raan,
                               eccentricity=eccentricity,
                               argperigee=argperigee,
                               meananomaly=meananomaly,
                               mm=mm)
    t1 = time.time()
    print(t1 - t0)
    t0 = time.time()
    ra, dec = config.tle.propagate_radec()
    t1 = time.time()
    print(t1 - t0)
    print('--')
    return ra, dec

def lnprior(theta, config):
    """
    Compute log prior - used to ensure walkers are exploring
    ranges defined by priors
    
    Parameters
    ----------
    theta : array-like
        Current set of parameters from MCMC
    config : mcmc.Config
        Config containing prior information
    
    Returns
    -------
    lnprior : float
        Log prior of current sample (0.0 | -np.inf)
    """
    for p in range(len(theta)):
        if theta[p] < config.priors['llim'][p] or theta[p] > config.priors['ulim'][p]:
            return -np.inf
    return 0.0

def lnlike(theta, config):
    """
    Compute log likelihood for proposed model

    Parameters
    ----------
    theta : array-like
        Current set of parameters from MCMC
    config : mcmc.Config
        Config containing prior and data information
    """
    # evaluate model
    ra_model, dec_model = model(theta, config)

    if True in np.isnan(ra_model) or True in np.isnan(dec_model):
        print('here :(')
        return -np.inf

    # compute log likelihood for RA
    lnlike_ra = (config.data['ra'] - ra_model) ** 2 / config.data['ra_sigma2'] + config.data['ra_lnsigma2']
    lnlike_ra = -0.5 * (np.sum(lnlike_ra) - np.log(len(config.data['ra']) + 1))

    # compute log likelihood for DEC
    lnlike_dec = (config.data['dec'] - dec_model) ** 2 / config.data['dec_sigma2'] + config.data['dec_lnsigma2']
    lnlike_dec = -0.5 * (np.sum(lnlike_dec) - np.log(len(config.data['dec']) + 1))

    # sum to return overall log likelihood
    return lnlike_ra + lnlike_dec

def lnprob(theta, config):
    """
    Log probability function - wraps lnprior and lnlike

    Parameters
    ----------
    theta : array-like
        Current set of parameters from MCMC
    config : mcmc.Config
        Config containing prior and data information
    """
    lp = lnprior(theta, config)
    print(lp)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, config)

def runSampler(config):
    """
    Set up sampler and run MCMC
    """
    sampler = emcee.EnsembleSampler(config.n_walkers,
                                    config.n_dim,
                                    lnprob,
                                    args=(config,))
    print('Running MCMC...')
    sampler.run_mcmc(config.init_pos,
                     config.n_steps,
                     rstate0=np.random.get_state())
    print('Done.')
    return sampler
