"""
General tasks for tlemcee MCMC algorithm
"""

import numpy as np
from datetime import datetime

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

    config.tle.modify_elements(mmdot=mmdot,
                               mmdot2=mmdot2,
                               drag=drag,
                               inclination=inclination,
                               raan=raan,
                               eccentricity=eccentricity,
                               argperigee=argperigee,
                               meananomaly=meananomaly,
                               mm=mm)

    return config.tle.propagate_radec()

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
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, config)

def getBestModel(sampler, burnin=None, save=False):
    """
    Obtain best model from a given sampler chain

    Parameters
    ----------
    sampler : mcmc.sampler.Sampler
        MCMC sampler containing chain information
    burnin : int, optional
        Burn-in phase for sampler chain (user-assessed)
        If given, use median method on sampler chain with burnin phase removed
        If None, use argmax method on entire sampler chain
        Default = None
    save : bool, optional
        Toggle to save best model TLE to file
        Default = False

    Returns
    -------
    best_model : st.tle.TLE
        Best model TLE
    """
    best_params = {}
    if burnin is None:
        # use argmax method
        best_pars_idx = np.unravel_index(np.argmax(sampler.sampler.lnprobability),
                                         (sampler.config.n_walkers,
                                          sampler.config.n_steps))
        best_pars = sampler.sampler.chain[best_pars_idx[0], best_pars_idx[1], :]

        for p, param in enumerate(sampler.config.priors['name']):
            best_params[param] = {'value': best_pars[p],
                                  'error': None}
    else:
        # use median method
        samples = sampler.sampler.chain[:, burnin:, :].reshape((-1,
                                                                sampler.config.n_dim))

        for p, param in enumerate(sampler.config.priors['name']):
            best_params[param] = {'value': np.median(samples[:, i]),
                                  'error': np.std(samples[:, i])}

    if 'mmdot' in best_params:
        mmdot = best_params['mmdot']['value']
    else:
        mmdot = None
    if 'mmdot2' in best_params:
        mmdot2 = best_params['mmdot2']['value']
    else:
        mmdot2 = None
    if 'drag' in best_params:
        drag = best_params['drag']['value']
    else:
        drag = None
    if 'inclination' in best_params:
        inclination = best_params['inclination']['value']
    else:
        inclination = None
    if 'raan' in best_params:
        raan = best_params['raan']['value']
    else:
        raan = None
    if 'eccentricity' in best_params:
        eccentricity = best_params['eccentricity']['value']
    else:
        eccentricity = None
    if 'argperigee' in best_params:
        argperigee = best_params['argperigee']['value']
    else:
        argperigee = None
    if 'meananomaly' in best_params:
        meananomaly = best_params['meananomaly']['value']
    else:
        meananomaly = None
    if 'mm' in best_params:
        mm = best_params['mm']['value']
    else:
        mm = None

    sampler.config.tle.modify_elements(mmdot=mmdot,
                                       mmdot2=mmdot2,
                                       drag=drag,
                                       inclination=inclination,
                                       raan=raan,
                                       eccentricity=eccentricity,
                                       argperigee=argperigee,
                                       meananomaly=meananomaly,
                                       mm=mm)
    if save:
        sampler.config.tle.save(sampler.out_dir)
    return sampler.config.tle
